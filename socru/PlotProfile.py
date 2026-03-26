"""
Circular genome plot generation for visualization.

This module creates circular pie chart visualizations of genome structure,
showing inter-operon fragments as wedges with their sizes proportional to
fragment length. Reversed fragments are indicated with hatching, and the
origin (dnaA) is labeled.

Classes:
    PlotProfile: Creates circular genome structure plots
"""

import matplotlib.pyplot as plt

from socru.SvgGenomePlot import generate_genome_svg, save_genome_svg

class PlotProfile:
    """
    Create circular visualization of genome structure.
    
    This class generates a pie chart representation of a bacterial genome
    showing inter-operon fragments arranged circularly. Each fragment is
    a wedge sized by its length. Reversed fragments get hatched patterns,
    and the origin of replication is labeled.
    
    Attributes:
        fragments (list): List of Fragment objects to plot
        output_file (str): Path where plot PDF will be saved
        verbose (bool): Enable verbose output
    """
    
    def __init__(self, fragments, output_file, verbose):
        """
        Initialize PlotProfile with fragments and output path.
        
        Args:
            fragments (list): Ordered list of Fragment objects
            output_file (str): Path to output PDF file
            verbose (bool): Enable verbose output
        """
        self.fragments = fragments
        self.output_file = output_file
        self.verbose = verbose
        
    def total_bases(self):
        """
        Calculate total genome size from all fragments.
        
        Returns:
            int: Sum of all fragment lengths
        """
        return sum([f.num_bases() for f in self.fragments])
    
    def create_plot(self):
        """
        Create and save circular genome structure plot.
        
        Generates a pie chart where:
        - Each wedge represents a fragment
        - Wedge size is proportional to fragment length
        - Reversed fragments have hatched pattern ('+')
        - Origin fragment is labeled with "- Ori"
        - White circle in center creates donut appearance
        
        Saves plot as PDF to output_file.
        """
        # Prepare fragment data for plotting
        names = [str(f.number) for f in self.fragments]
        size = [int(f.num_bases()) for f in self.fragments]
        reversed_fragments = [f.reversed_frag for f in self.fragments]
        dna_A = [f.dna_A for f in self.fragments]

        # Label the origin fragment
        for i in range(len(names)):
            if dna_A[i]:
                names[i] += " - Ori"

        # Create white circle for donut hole
        my_circle=plt.Circle( (0,0), 0.7, color='white')
        
        # Create pie chart with white borders between wedges
        piechart = plt.pie(size, labels=names, wedgeprops = { 'linewidth' : 7, 'edgecolor' : 'white' })
        
        # Add hatching pattern to reversed fragments
        for i in range(len(piechart[0])):
            if reversed_fragments[i]:
                piechart[0][i].set_hatch('+')
        
        # Add the donut hole and save
        p=plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig( self.output_file )

    def _build_fragment_dicts(self):
        """Convert Fragment objects to dicts for the SVG generator.

        Returns:
            list: List of fragment dictionaries.
        """
        frag_dicts = []
        for f in self.fragments:
            total_len = f.num_bases()
            coords = [(c[0], c[1]) for c in f.coords] if f.coords else [(0, total_len)]
            frag_dicts.append({
                "number": f.number,
                "reversed": f.reversed_frag,
                "length": total_len,
                "coords": coords,
                "is_dnaA": f.dna_A,
                "is_dif": f.dif,
            })
        return frag_dicts

    def _build_operon_dicts(self):
        """Infer operon positions from fragment boundaries.

        Operons sit between adjacent fragments. Their direction is derived
        from the operon_forward_end / operon_forward_start flags on the
        surrounding fragments.

        Returns:
            list: List of operon dictionaries.
        """
        operons = []
        for idx in range(len(self.fragments)):
            frag = self.fragments[idx]
            next_frag = self.fragments[(idx + 1) % len(self.fragments)]
            # The operon sits at the end of this fragment / start of next
            if frag.coords:
                end_coord = frag.coords[-1][1]
            else:
                end_coord = 0
            if next_frag.coords:
                start_coord = next_frag.coords[0][0]
            else:
                start_coord = end_coord
            direction = "forward" if frag.operon_forward_end else "reverse"
            operons.append({
                "start": end_coord,
                "end": start_coord,
                "direction": direction,
            })
        return operons

    def create_svg(self, svg_output_file, gs_type="", quality="GREEN", genome_name=""):
        """Create and save an SVG circular genome diagram.

        Args:
            svg_output_file (str): Path to output SVG file.
            gs_type (str): GS type label, e.g. "GS1.0".
            quality (str): Quality category: "GREEN", "AMBER", or "RED".
            genome_name (str): Genome filename or label.
        """
        frag_dicts = self._build_fragment_dicts()
        operon_dicts = self._build_operon_dicts()
        genome_length = self.total_bases()

        save_genome_svg(
            svg_output_file,
            fragments=frag_dicts,
            operons=operon_dicts,
            genome_length=genome_length,
            gs_type=gs_type,
            quality=quality,
            genome_name=genome_name,
        )