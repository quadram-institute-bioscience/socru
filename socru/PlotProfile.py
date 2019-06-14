import matplotlib.pyplot as plt

class PlotProfile:
    
    def __init__(self, fragments, output_file):
        self.fragments = fragments
        self.output_file = output_file
        
    def total_bases(self):
        return sum([f.num_bases() for f in self.fragments])
    
    # make this look pretty
    def create_plot(self):
        names = [str(f.number) for f in self.fragments]
        size = [int(f.num_bases()) for f in self.fragments]
        reversed_fragments = [f.reversed_frag for f in self.fragments]
        dna_A = [f.dna_A for f in self.fragments]
        
        for i in range(len(names)):
            if dna_A[i]:
                names[i] += " - Ori"

        my_circle=plt.Circle( (0,0), 0.7, color='white')
        
        piechart = plt.pie(size, labels=names, wedgeprops = { 'linewidth' : 7, 'edgecolor' : 'white' })
        
        # set hatching pattern
        for i in range(len(piechart[0])):
            if reversed_fragments[i]:
                piechart[0][i].set_hatch('+')
        p=plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig( self.output_file )
        
  