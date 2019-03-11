# classes for plotting
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from pymatgen.core.periodic_table import Element

mpl.rcParams["savefig.directory"] = os.chdir(os.getcwd())
mpl.rcParams["savefig.format"] = 'pdf'





#TODO: clean this class a bit
class PlotterPSE:

    def __init__(self, valuestoplot,counter_cations_env=None):
        """

        :param cationlist:
        :param valuestoplot:
        :param counter_cations_env: let us determine the number of environments for each cat
        """

        self.cationlist = valuestoplot.keys()
        self.valuestoplot = valuestoplot
        self.counter_cations_env=counter_cations_env


    def _get_fulfilled(self, okay, notokay):
        okayperc = np.float(okay) / np.float(okay + notokay)
        return okayperc

    def set_colormap(self,lowerlimit,upperlimit):
        self.fig.colorbar(self.img1, ax=self.ax).set_clim(lowerlimit, upperlimit)
        #self.ax.clim(lowerlimit, upperlimit)

    def get_plot(self, xlim=[1, 18], ylim=[1, 9],lowest_number_of_environments_considered=50):
        import matplotlib
        # TODO: cleanup
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        matplotlib.rcParams['font.size'] = 3
        mpl.rcParams['figure.dpi'] = 400
        matplotlib.rc("savefig", dpi=400)
        #font = {'family': 'normal',
        #        'size': 10}

        #matplotlib.rc('font', **font)

        PSE = np.full((int(ylim[1] - ylim[0] + 2),
                       int(xlim[1] - xlim[0] + 2)), np.nan)

        for cat in self.cationlist:
            if self.counter_cations_env is None:
                if not self.valuestoplot[cat][0]+self.valuestoplot[cat][1]<lowest_number_of_environments_considered:
                    if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                        PSE[Element(cat).row, Element(cat).group] = self._get_fulfilled(self.valuestoplot[cat][0],
                                                                                        self.valuestoplot[cat][1])
            else:
                if not self.counter_cations_env[cat]<lowest_number_of_environments_considered:
                    if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                        PSE[Element(cat).row, Element(cat).group] = self._get_fulfilled(self.valuestoplot[cat][0],
                                                                                        self.valuestoplot[cat][1])

        fig, ax = plt.subplots()
        mycm = plt.cm.get_cmap('cool', 100)
        img1 = ax.imshow(PSE, cmap=mycm)

        for cat in self.cationlist:
            if self.counter_cations_env is None:
                if not (self.valuestoplot[cat][0] + self.valuestoplot[cat][1] < lowest_number_of_environments_considered):
                    if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                        ax.text(Element(cat).group - 0.3, Element(cat).row + 0.15, cat)
            else:
                if not self.counter_cations_env[cat] < lowest_number_of_environments_considered:
                    if cat=='Np':
                        print(self.counter_cations_env[cat])
                    if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                        ax.text(Element(cat).group - 0.3, Element(cat).row + 0.15, cat)

        ax.set_ylim(np.float(ylim[1]) - 0.5, np.float(ylim[0] - 0.5))
        ax.set_xlim(np.float(xlim[0] - 0.5), np.float(xlim[1] - 0.5))
        ax.set_xticks(range(int(xlim[0]), int(xlim[1]), 1))
        ax.set_yticks(range(int(ylim[0]), int(ylim[1]), 1))
        current_cmap = mpl.cm.get_cmap()
        current_cmap.set_bad(color='white')
        self.fig=fig
        self.img1=img1
        self.ax=ax
        return plt
