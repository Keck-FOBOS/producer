"""
Simulate allocating FOBOS apertures to a randomly drawn set of targets.

Contains code originally written by 2021 Akamai intern, Brittany Ann
Ramos.

.. include:: ../include/links.rst
"""

from IPython import embed 

import numpy
from matplotlib import pyplot, patches, ticker

from . import scriptbase
from ..targets import random_targets
from ..plan import configure_observations


class AllocSim(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description='FOBOS Allocation Simulation', width=width)

        parser.add_argument('-d', '--density', nargs=3, type=float, default=[2.5, 40, 5],
                            help='Target density sampling: minimum, maximum, number of samples.  '
                                 'Density is sampled geometrically.')
        parser.add_argument('-s', '--sims', type=int, default=1,
                            help='Number of simulations to use for mean and standard deviation '
                                 'of trends.') 
        parser.add_argument('-m', '--mode', type=int, default=1,
                            help='The spectrograph mode.  All spectrographs are put in the same '
                                 'mode.   Use 1 for MOS mode, 2 for multi-IFU mode.')

        parser.add_argument('-e', '--embed', default=False, action='store_true',
                            help='Embed using IPython before completing the script.')

        return parser

    @staticmethod
    def main(args):

        rng = numpy.random.default_rng(99)

        ndens = int(args.density[2])
        nsim = args.sims
        density = numpy.geomspace(args.density[0], args.density[1], ndens)
        nobj = numpy.zeros((nsim, ndens), dtype=int)
        nobs = numpy.zeros((nsim, ndens), dtype=int)
        nalloc = numpy.empty((nsim, ndens), dtype=object)
        nap = numpy.empty((nsim, ndens), dtype=object)
        completeness = numpy.empty((nsim, ndens), dtype=object)
        efficiency = numpy.empty((nsim, ndens), dtype=object)
        mean_efficiency = numpy.empty((nsim, ndens), dtype=object)

        for j in range(nsim):
            for i in range(ndens):
                print(f'Density {i+1}/{ndens}', end='\r')

                objx, objy = random_targets(10., density=density[i], rng=rng)
                aptype = numpy.zeros(objx.size, dtype=int)

                n_in_fov, obs_obj, obs_nap, obs_ap, obs_mode \
                        = configure_observations(objx, objy, aptype)

                nobj[j,i] = n_in_fov
                nobs[j,i] = len(obs_obj)
                nalloc[j,i] = numpy.array([len(o) for o in obs_obj])
                nap[j,i] = numpy.array(obs_nap)
                completeness[j,i] = numpy.cumsum(nalloc[j,i])/nobj[j,i]
                efficiency[j,i] = nalloc[j,i]/nap[j,i]
                mean_efficiency[j,i] = numpy.cumsum(efficiency[j,i])/(numpy.arange(nobs[j,i])+1)

            print(f'Density {ndens}/{ndens}')

        obs_mask = numpy.empty((nsim, ndens), dtype=object)
        max_nobs = numpy.amax(nobs, axis=0)
        for j in range(nsim):
            for i in range(ndens):
                obs_mask[j,i] = numpy.zeros(max_nobs[i], dtype=bool)
                if nobs[j,i] < max_nobs[i]:
                    obs_mask[j,i][nobs[j,i]:] = True
                    nalloc[j,i] = numpy.append(nalloc[j,i], numpy.zeros(max_nobs[i] - nobs[j,i]))
                    nap[j,i] = numpy.append(nap[j,i], numpy.zeros(max_nobs[i] - nobs[j,i]))
                    completeness[j,i] = numpy.append(completeness[j,i],
                                                numpy.zeros(max_nobs[i] - nobs[j,i]))
                    efficiency[j,i] = numpy.append(efficiency[j,i], numpy.zeros(max_nobs[i] - nobs[j,i]))
                    mean_efficiency[j,i] = numpy.append(mean_efficiency[j,i],
                                                    numpy.zeros(max_nobs[i] - nobs[j,i]))

        m_completeness = numpy.empty(ndens, dtype=object)
        m_efficiency = numpy.empty(ndens, dtype=object)
        m_fom = numpy.empty(ndens, dtype=object)
    
        s_completeness = numpy.empty(ndens, dtype=object)
        s_efficiency = numpy.empty(ndens, dtype=object)
        s_fom = numpy.empty(ndens, dtype=object)
    
        for i in range(ndens):
            comp = numpy.ma.MaskedArray(numpy.stack(completeness[:,i]), mask=numpy.stack(obs_mask[:,i]))
            m_completeness[i] = numpy.ma.mean(comp, axis=0).filled(0.0)
            s_completeness[i] = numpy.ma.std(comp, axis=0).filled(0.0)

            arr = numpy.ma.MaskedArray(numpy.stack(efficiency[:,i]), mask=numpy.stack(obs_mask[:,i]))
            m_efficiency[i] = numpy.ma.mean(arr, axis=0).filled(0.0)
            s_efficiency[i] = numpy.ma.std(arr, axis=0).filled(0.0)

            arr = numpy.ma.sqrt(comp) \
                    * numpy.ma.MaskedArray(numpy.stack(mean_efficiency[:,i]), mask=numpy.stack(obs_mask[:,i]))
            m_fom[i] = numpy.ma.mean(arr, axis=0).filled(0.0)
            s_fom[i] = numpy.ma.std(arr, axis=0).filled(0.0)

        logformatter = ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(
                                        int(numpy.maximum(-numpy.log10(y),0)))).format(y))

        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))
        ax = fig.add_axes([0.2, 0.69, 0.6, 0.3])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.set_xlim(0.8, numpy.amax(nobs) + 0.2)
#        ax.set_ylim(0., 1.05)
        ax.set_ylim(0.01, 1.05)
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(logformatter)
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        for i in range(ndens):
#            ax.fill_between(numpy.arange(max_nobs[i])+1, m_completeness[i] + s_completeness[i],
#                            y2=m_completeness[i] - s_completeness[i],
#                            color=f'C{i}', lw=0, alpha=0.2, zorder=1)
            _x = numpy.arange(max_nobs[i])+1
            _y = 1 - m_completeness[i]
            ax.scatter(_x, _y, color=f'C{i}', marker='.',  s=30, lw=0, zorder=3)
            ax.plot(_x, _y, color=f'C{i}', zorder=3, label=r'$\rho$' + f' = {density[i]:.1f}')
        ax.text(-0.1, 0.5, '1-Completeness', ha='center', va='center', transform=ax.transAxes,
                rotation='vertical')
        ax.legend()

        ax = fig.add_axes([0.2, 0.38, 0.6, 0.3])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.set_xlim(0.8, numpy.amax(nobs) + 0.2)
#        ax.set_ylim(0., 1.05)
        ax.set_ylim(0.01, 1.05)
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(logformatter)
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        for i in range(ndens):
#            ax.fill_between(numpy.arange(max_nobs[i])+1, m_efficiency[i] + s_efficiency[i],
#                            y2=m_efficiency[i] - s_efficiency[i],
#                            color=f'C{i}', lw=0, alpha=0.2, zorder=1)
            _x = numpy.arange(max_nobs[i])+1
            _y = m_efficiency[i]
            ax.scatter(_x, _y, color=f'C{i}', marker='.', s=30, lw=0, zorder=3)
            ax.plot(_x, _y, color=f'C{i}', zorder=3)
        ax.text(-0.1, 0.5, 'Efficiency', ha='center', va='center', transform=ax.transAxes,
                rotation='vertical')

        ax = fig.add_axes([0.2, 0.07, 0.6, 0.3])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.set_xlim(0.8, numpy.amax(nobs) + 0.2)
#        ax.set_ylim(0., 1.05)
        ax.set_ylim(0.1, 1.05)
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(logformatter)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        for i in range(ndens):
#            ax.fill_between(numpy.arange(max_nobs[i])+1, m_fom[i] + s_fom[i],
#                            y2=m_fom[i] - s_fom[i],
#                            color=f'C{i}', lw=0, alpha=0.2, zorder=1)
            _x = numpy.arange(max_nobs[i])+1
            _y = m_fom[i]
            ax.scatter(_x, _y, color=f'C{i}', marker='.', s=30, lw=0, zorder=3)
            ax.plot(_x, _y, color=f'C{i}', zorder=3)
        ax.text(-0.1, 0.5, 'Figure-of-Merit', ha='center', va='center', transform=ax.transAxes,
                rotation='vertical')
        ax.text(0.5, -0.15, 'Pointing Number', ha='center', va='center', transform=ax.transAxes)

        ofile = None
        ofile = 'example_allocsim.png'
        if ofile is None:
            pyplot.show()
        else:
            fig.canvas.print_figure(ofile, bbox_inches='tight')
        fig.clear()
        pyplot.close(fig)













