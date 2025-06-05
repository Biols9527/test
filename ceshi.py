#!/usr/bin/env python3

import argparse
import csv
import networkx as nx
from pathlib import Path
from typing import Dict, Set, List, Optional, Any # Added Any for pickle load
import textwrap
import synphoni.utils as su
import synphoni.graph_analysis as sg
from synphoni.logo import logo_ASCII
import pickle
from dataclasses import dataclass
import sys # For sys.getsizeof
import gc # For garbage collection

# --- memory_profiler import and decorator ---
# Ensure you have memory_profiler installed: pip install memory_profiler
# To use it, run: python -m memory_profiler your_script_name.py [your_args]
try:
    from memory_profiler import profile
except ImportError:
    # Define a dummy profile decorator if memory_profiler is not installed
    # This allows the script to run without memory_profiler, just without profiling
    def profile(func):
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper
    print("Memory profiler not found. Install with 'pip install memory_profiler' for detailed memory usage.")

@dataclass
class SynphoniConfig:
    """Configuration class to store all parameters"""
    min_len: int
    min_shared: float
    clique_size: int
    min_community_coverage: float
    chrom_clustering_method: str

class SynphoniProcessor:
    """Main class to handle SYNPHONI processing"""

    @profile
    def __init__(self, config: SynphoniConfig):
        print(f"[{self.__class__.__name__}] Initializing with config: {config}")
        self.config = config
        self.chrom_dict: Dict[Any, Any] = {}
        self.G_og: Optional[nx.Graph] = None # Type hint for NetworkX graph
        self.species_ls: List[Any] = []
        self.ortho: Dict[Any, Any] = {}
        print(f"[{self.__class__.__name__}] Initialized.")

    @profile
    def load_data(self, filtered_graph_path: Path, chrom_data_path: Path) -> None:
        """Load required data files"""
        print(f"[{self.__class__.__name__}] Starting data loading...")

        # 1. Load filtered_graph
        print(f"  Loading filtered graph from: {filtered_graph_path}")
        try:
            with open(filtered_graph_path, "rb") as fhin:
                self.G_og = pickle.load(fhin)
            if self.G_og is not None:
                print(f"  SUCCESS: Loaded G_og. Type: {type(self.G_og)}.")
                if isinstance(self.G_og, nx.Graph):
                    print(f"    G_og info: Nodes={self.G_og.number_of_nodes()}, Edges={self.G_og.number_of_edges()}")
                else:
                    print(f"    G_og memory size (approx): {sys.getsizeof(self.G_og)} bytes")
            else:
                print(f"  WARNING: G_og loaded as None from {filtered_graph_path}")

        except (FileNotFoundError, pickle.UnpicklingError, EOFError) as e: # Added EOFError
            raise RuntimeError(f"Error loading filtered graph '{filtered_graph_path}': {str(e)}")
        print("-" * 30)

        # 2. Load chrom_data
        print(f"  Loading chromosome data from: {chrom_data_path}")
        try:
            self.chrom_dict = su.load_chrom_data(filepath=str(chrom_data_path)) # Ensure filepath is str if lib expects
            print(f"  SUCCESS: Loaded chrom_dict. Type: {type(self.chrom_dict)}. Number of OGs: {len(self.chrom_dict)}")
            # print(f"    chrom_dict memory size (approx): {sys.getsizeof(self.chrom_dict)} bytes") # sys.getsizeof might be misleading for complex dicts
        except (FileNotFoundError, pickle.UnpicklingError, EOFError, Exception) as e: # Catch generic Exception from su.load_chrom_data
            raise RuntimeError(f"Error loading chromosome data '{chrom_data_path}' using su.load_chrom_data: {str(e)}")
        print("-" * 30)

        # 3. Build ortho dictionary
        print("  Building ortho dictionary...")
        if not self.chrom_dict:
            print("  WARNING: chrom_dict is empty. Ortho dictionary will be empty.")
            self.ortho = {}
        else:
            try:
                # This comprehension can be memory intensive if chrom_dict is huge.
                # Monitor memory usage here carefully.
                self.ortho = {
                    acc: og_key
                    for og_key, og_data in self.chrom_dict.items()
                    for species, species_data in og_data.items()
                    for chromo, chromo_data in species_data.items()
                    for acc in chromo_data.keys() # Assuming chromo_data is a dict with acc as keys
                }
                print(f"  SUCCESS: Built ortho dictionary. Number of accessions mapped: {len(self.ortho)}")
            except Exception as e:
                raise RuntimeError(f"Error building ortho dictionary: {str(e)}")
        print("-" * 30)

        # 4. Build species_ls
        print("  Building species list...")
        if not self.chrom_dict:
             print("  WARNING: chrom_dict is empty. Species list will be empty.")
             self.species_ls = []
        else:
            try:
                self.species_ls = list(set(su.flatten([list(og_data.keys()) for og_key, og_data in self.chrom_dict.items()])))
                # Adjusted based on typical structure: chrom_dict[og] gives dict of species.
                # If su.flatten expects list of lists of species:
                # self.species_ls = list(set(su.flatten([[species for species in og_data.keys()] for og_data in self.chrom_dict.values()])))
                print(f"  SUCCESS: Built species list. Number of unique species: {len(self.species_ls)}")
            except Exception as e:
                raise RuntimeError(f"Error building species list: {str(e)}")

        print(f"[{self.__class__.__name__}] Data loading finished.")
        print("=" * 50)


    @profile
    def process_communities(self, og_communities_path: Path, output_prefix: str) -> Dict:
        """Process OG communities and write results"""
        print(f"[{self.__class__.__name__}] Starting community processing...")
        block_ids: Dict[Any, Any] = {}

        # 1. Prepare output paths
        synt_path_str = f"{output_prefix}.len{self.config.min_len}.ol{self.config.min_shared}.synt"
        multi_sp_path_str = f"{output_prefix}.len{self.config.min_len}.ol{self.config.min_shared}.clusters"
        print(f"  Synteny output will be written to: {synt_path_str}")
        print(f"  Clusters output will be written to: {multi_sp_path_str}")
        print("-" * 30)

        # 2. Read communities
        # This part can consume a lot of memory if the CSV is large, as it reads all communities into a list.
        print(f"  Reading OG communities from: {og_communities_path}")
        og_commus: List[Set[str]] = []
        try:
            with open(og_communities_path, "r") as f:
                reader = csv.reader(f)
                # Consider processing row by row if memory is an issue and sorting all at once is the bottleneck
                # For now, keeping original logic:
                og_commus_temp = [set(row) for row in reader]
            print(f"  SUCCESS: Read {len(og_commus_temp)} communities. Now sorting...")
            og_commus = sorted(og_commus_temp, key=len, reverse=True)
            del og_commus_temp # Free memory from the unsorted list
            gc.collect() # Explicitly ask for garbage collection
            print(f"  SUCCESS: Sorted {len(og_commus)} communities by length (descending).")
        except FileNotFoundError:
            raise RuntimeError(f"OG communities file not found: {og_communities_path}")
        except Exception as e:
            raise RuntimeError(f"Error reading or sorting OG communities: {str(e)}")
        print("-" * 30)

        # 3. Process each community
        print(f"  Starting to process {len(og_commus)} communities...")
        total_communities = len(og_commus)
        if total_communities == 0:
            print("  No communities to process.")
            return block_ids

        try:
            with open(synt_path_str, "w", newline='') as synt_h, open(multi_sp_path_str, "w", newline='') as multi_sp_h:
                synt_w = csv.writer(synt_h, delimiter="\t")
                m_sp_w = csv.writer(multi_sp_h, delimiter="\t")

                for i, current_commu in enumerate(og_commus):
                    if (i + 1) % 10 == 0 or i == 0 or i == total_communities -1 : # Print progress
                        print(f"    Processing community {i+1}/{total_communities} (size: {len(current_commu)} OGs)...")

                    # Delegate to _process_single_community
                    # The update to block_ids happens here. If block_ids grows very large, it's a concern.
                    new_ids = self._process_single_community(
                        current_commu, synt_w, m_sp_w, block_ids # Pass current block_ids for context if needed by write_blocks
                    )
                    if new_ids: # Assuming write_blocks returns only NEWLY created block IDs
                        block_ids.update(new_ids)

                    # Optional: Monitor block_ids size periodically
                    if (i + 1) % 100 == 0:
                        print(f"      Current block_ids count: {len(block_ids)}")
                        # print(f"      block_ids memory (approx): {sys.getsizeof(block_ids)} bytes") # Misleading for complex dicts
                        # gc.collect() # Optional: collect garbage more frequently if memory is tight

        except IOError as e:
            raise RuntimeError(f"Error opening or writing to output files: {str(e)}")
        except Exception as e: # Catch other errors during processing loop
            raise RuntimeError(f"Error during processing of community {i+1 if 'i' in locals() else 'unknown'}: {str(e)}")

        print(f"[{self.__class__.__name__}] Community processing finished. Total block IDs generated: {len(block_ids)}")
        print("=" * 50)
        return block_ids

    @profile
    def _process_single_community(self,
                                community: Set[str],
                                synt_writer: csv.writer,
                                multi_sp_writer: csv.writer,
                                existing_block_ids: Dict # Pass this for context if sg.write_blocks needs it
                               ) -> Dict:
        """Process a single OG community. Returns NEWLY created block IDs."""
        # print(f"      _process_single_community: Processing community with {len(community)} OGs.")

        # 1. Get scaffold locations
        current_commu_scaffolds: Optional[Dict] = None # Initialize
        try:
            current_commu_scaffolds = sg.genome_location_ogs(
                og_community=community,
                chrom_data=self.chrom_dict,
                species_list=self.species_ls,
                orthology=self.ortho,
                min_og_commu=self.config.min_community_coverage
            )
        except Exception as e:
            print(f"      ERROR in sg.genome_location_ogs for community {list(community)[:5]}...: {str(e)}")
            return {} # Return empty dict on error

        if not current_commu_scaffolds:
            # print("        No significant scaffold locations found for this community.")
            return {}

        # 2. Build protoblock graph
        protoblock_graph: Optional[nx.Graph] = None # Initialize
        try:
            protoblock_graph = sg.og_info_to_graph(
                genome_location_orthogroups=current_commu_scaffolds,
                fullgraph_ogs_filt=self.G_og,
                min_len=self.config.min_len,
                min_shared=self.config.min_shared
            )
        except Exception as e:
            print(f"      ERROR in sg.og_info_to_graph for community {list(community)[:5]}...: {str(e)}")
            return {} # Return empty dict on error

        # 3. Write blocks if graph exists
        if protoblock_graph is not None and protoblock_graph.number_of_nodes() > 0:
            try:
                # Assuming sg.write_blocks returns a dictionary of *newly created* block IDs for this community
                # or modifies existing_block_ids and returns the relevant subset or all of it.
                # For clarity, it's better if it returns only what's new or relevant to this call.
                newly_created_ids = sg.write_blocks(
                    blocks_writer=synt_writer,
                    multi_sp_writer=multi_sp_writer,
                    genome_location_ogs_dict=current_commu_scaffolds,
                    og_info_graph=protoblock_graph,
                    k_perco=self.config.clique_size,
                    known_dict=existing_block_ids, # Pass the main dict if the function needs to check for existing IDs
                    method=self.config.chrom_clustering_method
                )
                # If sg.write_blocks modifies existing_block_ids in place and returns None or a status,
                # then this part needs adjustment. The current code assumes it returns the new IDs.
                return newly_created_ids if newly_created_ids is not None else {}
            except Exception as e:
                print(f"      ERROR in sg.write_blocks for community {list(community)[:5]}...: {str(e)}")
                return {}
        else:
            # print("        Protoblock graph is None or empty. No blocks to write.")
            return {}
        # No explicit del needed here for current_commu_scaffolds, protoblock_graph
        # as they are local to this function and will be garbage collected when they go out of scope.
        # However, if memory profiling shows they are very large and persist, explicit del could be tried.


def main():
    print("Starting SYNPHONI Step 4 script...")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(f"""\
            {logo_ASCII()}
            Step 4 of the SYNPHONI pipeline:
            Extract syntenic blocks from orthogroup communities.
            This script processes orthogroup communities to identify syntenic blocks
            across multiple species based on phylogeny and ortholog network inference.
            Run with 'python -m memory_profiler {sys.argv[0]} [args...]' for memory profiling.
            """)
    )
    # ... (rest of your argparse setup is good) ...
    parser.add_argument("og_communities", help="CSV file containing orthogroup communities (output from step 3)", type=Path)
    parser.add_argument("-g", "--filtered_graph", help="Filtered graph in gpickle format (output from step 3)", type=Path, required=True)
    parser.add_argument("-c", "--chrom_data", help="Pickle file containing chromosome data (generated in step 1)", type=Path, required=True)
    parser.add_argument("-o", "--output", help="Prefix for output files (.synt and .clusters will be added)", type=str, required=True)
    parser.add_argument("-l", "--min_len", help="Minimum number of ortholog occurrences required", default=3, type=int)
    parser.add_argument("-s", "--min_shared", help="Minimum overlap coefficient between scaffolds (0.0-1.0)", default=0.5, type=float)
    parser.add_argument("-k", "--clique_size", help="Minimum size of multi-species blocks to retain", default=3, type=int)
    parser.add_argument("-r", "--min_community_coverage", help="Minimum percentage of orthogroups required (0.0-1.0)", default=0.3, type=float)
    parser.add_argument("-m", "--chrom_clustering_method", help="Method for clustering chromosomes", default="k_clique", choices={"k_clique", "leiden"}, type=str)

    args = parser.parse_args()
    print(f"Arguments parsed: {args}")
    print("=" * 50)

    # Validate arguments
    print("Validating input file paths...")
    if not args.og_communities.exists():
        print(f"FATAL: Communities file not found: {args.og_communities}")
        sys.exit(1)
    if not args.filtered_graph.exists():
        print(f"FATAL: Filtered graph file not found: {args.filtered_graph}")
        sys.exit(1)
    if not args.chrom_data.exists():
        print(f"FATAL: Chromosome data file not found: {args.chrom_data}")
        sys.exit(1)
    print("SUCCESS: All input files found.")
    print("=" * 50)

    # Create configuration
    config = SynphoniConfig(
        min_len=args.min_len,
        min_shared=args.min_shared,
        clique_size=args.clique_size,
        min_community_coverage=args.min_community_coverage,
        chrom_clustering_method=args.chrom_clustering_method
    )

    # Initialize and run processor
    print("Initializing SynphoniProcessor...")
    try:
        processor = SynphoniProcessor(config)
        print("-" * 30)

        print("Calling processor.load_data()...")
        processor.load_data(args.filtered_graph, args.chrom_data)
        print("-" * 30)

        print("Calling processor.process_communities()...")
        final_block_ids = processor.process_communities(args.og_communities, args.output)
        print("-" * 30)

        print(f"Script finished successfully. Total unique block IDs in final collection: {len(final_block_ids)}")

    except RuntimeError as e: # Catch custom runtime errors from processor
        print(f"RUNTIME ERROR during processing: {str(e)}")
        sys.exit(1)
    except Exception as e:
        print(f"UNEXPECTED ERROR during processing: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
