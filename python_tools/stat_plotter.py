
import matplotlib
matplotlib.use('Agg')
import psycopg2
import math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import ac_significance as ac_sig
from fcst_stat_database_access import experiment as fcst_db
import sys
var2Name = {}
var2Name['h'] = 'Height'
var2Name['t'] = 'Temperature'
var2Name['q'] = 'Specific Humidity'
var2Name['u'] = 'U Wind'
var2Name['v'] = 'V Wind'
dom2Long = {}
dom2Long['s.hem'] = 'Southern Hemisphere'
dom2Long['n.hem'] = 'Northern Hemisphere'
dom2Long['tropics'] = 'Tropics'
dom2Long['global'] = 'Global'
stat2Long ={}
stat2Long['cor'] = 'Anomaly Correlation'
stat2Long['rms'] = 'Root Mean Squared Error'


class plotter_2d:
    
   
    def plotter_2d(self,\
                    startdate=2018080100,\
                    enddate=2018093000,\
                    expid1='bmk_o4_nf.21z',\
                    expid2='bmk_ctl_joff.21z',\
                    expid1_name='9.6 $\mu$m',\
                    expid2_name='control',\
                    variables=['h','t','u','v','q'],\
                    domains=['n.hem','s.hem','tropics'],\
                    sigplot=True,\
                    verbose=False,\
                    straightmean=False,\
                    confidence=0.95,\
                    title=None,\
                    verif1='gmao',\
                    verif2='gmao',\
                    dblevs=None,\
                    maxlev=1000.,\
                    minlev=100.,\
                    levs=None,\
                    zero_t_zero=True,\
                    statistic='cor',\
                    graphicOutput='.png'):
        """
        Main Plotter which will loop through various variables and domains. For each variable and domain
        it will produce a 2D (3D if you count color) pcolor plot of the difference between the forecasts (expid1 minus expid2).
        """
        if levs is None:
            levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])
        else:
            levs = np.array(levs)
        self.confidence = confidence
        self.experiment1 = fcst_db(expid1, startdate=startdate, enddate=enddate) 
        self.experiment2 = fcst_db(expid2, startdate=startdate, enddate=enddate) 
        self.straightmean = straightmean 
        self.verbose = verbose
        self.zero_t_zero = zero_t_zero
        self.expid1_name = expid1_name
        self.expid2_name = expid2_name 
        self.sigplot = sigplot      
        self.ac = ac_sig.ac_significance(confidence=confidence)
        self.graphicOutput = graphicOutput 
        if (dblevs): levs = experiment1.get_pressure_levels()
        self.levs = levs

        print(levs.size,'Levels=',levs)
        steps = np.array([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120])
        steps = np.array([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120])
       
        pltlevs = (levs + np.roll(levs,1))/2.
        pltlevs[0] = levs[0] - (pltlevs[1]-levs[0])
        pltlevs = np.append(pltlevs,levs[-1] - (levs[-2]-levs[-1])*0.5)
        pltsteps = np.array([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132]) - 6

        (levs2d, stps2d) = self.init_2d_arrays(steps, levs, pltlevs, pltsteps)
        self.steps = steps
        self.levs = levs 
        self.levs2d = levs2d
        self.stps2d = stps2d
        for var in variables:
            for dom in domains:
                self.populate_stats(dom, statistic, var)
                self.draw_pcolor(pltsteps, pltlevs, statistic, var, dom,\
                                 self.experiment1.verify, self.experiment2.verify)


    def init_2d_arrays(self, steps, levs, pltlevs, pltsteps):
        """
        Helper function which initializes arrays with zeros, and populates model levels and time steps.
        """
        self.corarr = np.zeros((steps.size,levs.size))
        self.pctarr = np.zeros((steps.size,levs.size))
        self.sigarr = np.zeros((steps.size,levs.size))
        self.exp1arr = np.zeros((steps.size,levs.size))
        self.exp2arr = np.zeros((steps.size,levs.size))
        levs2d=np.zeros((pltsteps.size,pltlevs.size))
        stps2d=np.zeros((pltsteps.size,pltlevs.size))

        levs2d[:,0:pltlevs.size] = pltlevs[0:pltlevs.size].reshape(1,pltlevs.size)    
        stps2d[0:pltsteps.size,:] = pltsteps[0:pltsteps.size].reshape(pltsteps.size,1)
      

        return (levs2d, stps2d) 

    def populate_stats(self, domain, statistic, variable):
        """
        Helepr fuction which populates the plot arrays for a given domain, statistic and variable. 
        """
        print('Populating Stats for {} {} {}'.format(statistic, dom2Long[domain], var2Name[variable]))
        for i in np.arange(self.steps.size):
            for j in np.arange(self.levs.size):
                cstp = self.steps[i]
                clev = self.levs[j]
                exp1 = self.experiment1.get_stats(clev, cstp, domain, statistic, variable)             
                exp2 = self.experiment2.get_stats(clev, cstp, domain, statistic, variable)
                n = exp2.size
                expid1 = self.experiment1.experiment
                expid2 = self.experiment2.experiment
                startdate = self.experiment1.startdate
                enddate = self.experiment1.enddate
                if(self.verbose): print('Starting...',cstp,clev,domain,statistic,expid1,expid2,variable,startdate,enddate,'size=',n,'confidence=',self.confidence)

                if(statistic=='cor'):ztdiff, sigrange, sig, expcor1, expcor2 = self.ac.get_ztran_diff(exp1, exp2)
                else: ztdiff, sigrange, sig, expcor1, expcor2 = self.ac.get_tran_diff(exp1, exp2)
                if (not self.straightmean):
                    self.exp1arr[i,j] = expcor1
                    self.exp2arr[i,j] = expcor2
                else:
                    self.exp1arr[i,j] = np.mean(exp1)
                    self.exp2arr[i,j] = np.mean(exp2)


                if (self.verbose): print('  ztr Mean Diff, ztr range:',ztdiff,sigrange)
                self.corarr[i,j] = ztdiff
                if(abs(ztdiff)>0.0): 
                    self.sigarr[i,j] = False if ztdiff < sigrange[0] and ztdiff > sigrange[1] else True
                else:
                    self.sigarr[i,j] = False
                self.cnt = exp2.size

    def draw_pcolor(self, pltsteps, pltlevs, statistic, var, dom, verif1, verif2):
        """
        Helper function which actually does the plotting.  
        """
        fig= plt.figure(figsize=(6,5)) 
        ax = fig.add_subplot(111)

        plt.subplots_adjust(top=0.85, left=0.18)
        if (self.zero_t_zero): self.corarr[0,:] = 0.0
        mx = np.max(self.corarr)
        mn = np.min(self.corarr)
        rng = max( ( np.abs(mx), np.abs(mn) ) )
       
        if(statistic == 'cor'):
            cm = 'PiYG'
        else:
            cm = 'PiYG_r'
        fg = plt.pcolor(self.stps2d, self.levs2d, self.corarr, cmap = cm, vmin = rng*-1.0, vmax = rng)
        
        plt.xlim([0 - 6,120 + 6])
        plt.ylim([pltlevs[0],pltlevs[-1]])           

        plt.yscale('log')
            
        plt.yticks(self.levs, self.levs)
        plt.xticks(self.steps)
            
        plt.xlabel('Forecast Hour')
        plt.ylabel('Pressure (hPa)')
        startdate = self.experiment1.startdate
        enddate = self.experiment1.enddate
        expid1 = self.experiment1.experiment
        expid2 = self.experiment2.experiment 
        n = self.cnt
        #make strings to make date look prettier.
        sd = '{}'.format(self.experiment1.startdate)
        ed = '{}'.format(self.experiment1.enddate)
        startdate_fmt = sd[0:4]+'.'+sd[4:6]+'.'+sd[6:8]+' '+sd[8:10]+' UTC' 
        enddate_fmt = ed[0:4]+'.'+ed[4:6]+'.'+ed[6:8]+' '+ed[8:10]+' UTC'
        startdate,enddate = sd,ed
 

        if (statistic =='cor'):
            outtitle = "Anomaly Correlation {} Difference {}\n{} minus {}\n{} - {}, {} Forecasts".format(var2Name[var], dom2Long[dom],self.expid1_name, self.expid2_name, startdate_fmt, enddate_fmt, n)
        elif(statistic =='rms'):
            outtitle = "RMSE {} Difference {}\n{} minus {}\n{} - {}, {} Forecasts".format(var2Name[var], dom2Long[dom], self.expid1_name, self.expid2_name, startdate_fmt, enddate_fmt, dom2Long[dom], n) 
        plt.title(outtitle, fontsize=10)
        if (self.sigplot):
            self.ac.oplot_sig_hatch(ax, self.sigarr, pltsteps, pltlevs)
        cbar= plt.colorbar()
        cbar.set_label('Difference')
        fn = '{}_z.{}.{}-{}.{}-{}.{}.verif-{}.{}'.format(statistic, var, expid1, expid2, startdate, enddate, dom, verif1, verif2)
        print('Writing: ',fn + self.graphicOutput) 
        fig.savefig(fn + self.graphicOutput)
        plt.close(fig)


class plotter_1d:
    def plotter_1d(self,\
                    startdate=2018080100,\
                    enddate=2018093000,\
                    exp_ids=['bmk_o4_nf.21z'],\
                    control='bmk_ctl_joff.21z',\
                    exp_names=['9.6 $\mu$m'],\
                    control_name='Control',\
                    variables=['h','t','u','v','q'],\
                    domains=['n.hem','s.hem','tropics'],\
                    verbose=False,\
                    confidence=0.95,\
                    verif1='gmao',\
                    verif2='gmao',\
                    dblevs=None,\
                    maxlev=1000.,\
                    minlev=100.,\
                    levs=None,\
                    statistic='cor',\
                    straightmean=False,\
                    ncol = 3,\
                    graphicOutput='.png'):
        """
        Main Plotter for 1d which will loop through various variables and domains. For each variable and domain
        """
        if levs is None:
            levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])
        else:
            levs = np.array(levs)
        self.confidence = confidence
        self.experiments = {}
        self.experiment_names = {}
        self.straightmean = straightmean
        self.control = fcst_db(control, startdate=startdate, enddate=enddate)
        for e in exp_ids:
            self.experiments[e] = fcst_db(e, startdate=startdate, enddate=enddate)
        for i,e in enumerate(exp_ids):
            self.experiment_names[e] = exp_names[i]
        self.exp_ids = exp_ids
        self.control_name = control_name 
        self.verbose = verbose
        self.ac = ac_sig.ac_significance(confidence=confidence)

        self.graphicOutput = graphicOutput 
        if (dblevs): levs = experiment1.get_pressure_levels()
        self.levs = levs

        print(levs.size,'Levels=',levs)
        steps = np.array([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120])
        steps = np.array([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120])

        self.init_arrays(steps,exp_ids)
        self.steps = steps
        self.levs = levs 
        for var in variables:
            for dom in domains:
                for l in levs:
                    self.populate_stats(exp_ids, dom, statistic, var,l)
                    self.draw_plot(l, statistic, var, dom, self.control.verify, self.control.verify,ncol=ncol)
    def init_arrays(self, steps, exp_ids):
        """
        Helper function which initializes arrays with zeros, and populates model levels and time steps.
        """
        self.corarr = {}
        self.pctarr = {}
        self.sigarr = {}
        self.exparr = {}
        for e in exp_ids:
            self.corarr[e] = np.zeros([steps.size])
            self.pctarr[e] = np.zeros([steps.size])
            self.sigarr[e] = np.zeros([steps.size,2])
            self.exparr[e] = np.zeros([steps.size])
        self.controlarr = np.zeros(steps.size)

    def populate_stats(self, experiment_keys, domain, statistic, variable, lev):
        """
        Helepr fuction which populates the plot arrays for a given domain, statistic and variable. 
        """
        print('Populating Stats for {} {} {}'.format(statistic, dom2Long[domain], var2Name[variable]))
        exp = {}
        for i in np.arange(self.steps.size):
            cstp = self.steps[i]
            clev = lev
            for e in experiment_keys:
                exp[e] = self.experiments[e].get_stats(clev, cstp, domain, statistic, variable) 
            control = self.control.get_stats(clev, cstp, domain, statistic, variable)
            n = control.size
            expid_control = self.control.experiment
            startdate = self.control.startdate
            enddate = self.control.enddate
            if(self.verbose): print('Starting...',cstp,clev,domain,statistic,expid_control,variable,startdate,enddate,'size=',n,'confidence=',self.confidence)
            for e in experiment_keys:
                if(statistic=='cor'): ztdiff, sigrange, sig, expcor1, expcor2 = self.ac.get_ztran_diff(exp[e], control)
                else: ztdiff, sigrange, sig, expcor1, expcor2 = self.ac.get_tran_diff(exp[e], control)
                if (not self.straightmean):
                    self.exparr[e][i] = expcor1
                    self.controlarr[i] = expcor2
                else:
                    self.exparr[e][i]  = np.mean(exp[e])
                    self.controlarr[i] = np.mean(control)


                if (self.verbose): print('  ztr Mean Diff, ztr range:',ztdiff,sigrange)
                self.corarr[e][i] = ztdiff
                self.sigarr[e][i,0] = sigrange[0]
                self.sigarr[e][i,1] = sigrange[1] 
            self.cnt = control.size

    def draw_plot(self, lev, statistic, var, dom, verif1, verif2, boxes_or_bars='boxes',ncol=3):
        """
        Helper function which actually does the plotting.  
        """
        #alright 15 colors, if you're comparing 15 experiments or more, you're probably smoking something.
        color_list = ['green','red','blue','orange','purple',\
                      'brown','pink','gray','olive','cyan',\
                      'gold','wheat','peachpuff','khaki','fuchsia']

        # Here we go! start with 2 subplots the statistic (top) difference (bottom)
        fig, axs = plt.subplots(2,1)

        #book keeping for colors and tracking plot handles to use fig.legend instead of ax.legend.
        experimentColors = {}
        pltHandles = []
        
        #start with top plot which is the overlay of statistics plotted
        pltHandles.append(axs[0].plot(self.steps,self.controlarr,color='black',label=self.control_name,marker='o',linewidth=0.5,markersize=2.0)[0])
        experimentColors['control'] = 'black'

        # Add Plot N  experiments provided by user. 
        for i,e in enumerate(self.exp_ids):
            pltHandles.append( axs[0].plot(self.steps,self.exparr[e],label=self.experiment_names[e],marker='o',color=color_list[i],linewidth=0.5,markersize=2.0)[0])
            experimentColors[e] = color_list[i]            

        # set ticks and drop xtick labels for top plot
        axs[0].set_xticks(self.steps)
        axs[0].set_xticklabels(" ")

        #Label the y-axis with long name for statistic (e.g., Anomaly Correlation for cor)
        axs[0].set_ylabel(stat2Long[statistic])

        # start "bonus legend" where we print out the statistic for the last day. 
        #sort which experiment with largest end value 
        endVals = {}
        endVals['control'] = self.controlarr[-1]
        for e in self.exp_ids: endVals[e] = self.exparr[e][-1]
        evSorted = {k: v for k, v in sorted(endVals.items(), key=lambda item: item[1])}
        offset=-1
        bottomVal = list(evSorted.keys())[0]
        for e in list(evSorted.keys()):
            if(endVals[e]>0.01):
                axs[0].annotate('{:10.4f}'.format(endVals[e]), (self.steps[-1],endVals[bottomVal]),xytext=(12,offset),textcoords="offset points", fontsize=8,color=experimentColors[e])
            else:
                axs[0].annotate('{:.4e}'.format(endVals[e]), (self.steps[-1],endVals[bottomVal]),xytext=(17,offset),textcoords="offset points", fontsize=8,color=experimentColors[e])
            offset+=8
        # end bonus end point value "legend"

        # start difference plot bit
        for i,e in enumerate(self.exp_ids):
            axs[1].plot(self.steps,self.corarr[e],label=self.experiment_names[e],marker='o',color=color_list[i],linewidth=0.5, markersize=2.0)
            if(boxes_or_bars == 'bars'): axs[1].errorbar(self.steps,np.zeros(self.steps.size),yerr=np.abs(self.sigarr[e]).T,fmt='none',elinewidth=0,capsize=10,color=color_list[i]) 
            elif(boxes_or_bars =='boxes'): self.make_significance_boxes(axs[1], self.steps, np.zeros(self.steps.size), 0.4*(self.steps[1]-self.steps[0])*np.ones([2,self.steps.size]), np.abs(self.sigarr[e][:,:]).T,edgecolor=color_list[i] )
            else: sys.exit('unknown selection for boxes or bars.')
        axs[1].set_xticks(self.steps)
        axs[1].set_ylabel('Difference')
        axs[1].set_xlabel('Forecast Hour')
        # make sure we don't crowd out the y-axis with long decimal places force to use scientific notation past 3 decimal places.
        axs[1].ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
        axs[0].ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
        # Draw and x-axis for the difference plot
        axs[1].axhline(color='black')  

        # funk little piece to give enough space for a legend at the bottom.
        # variable size depending how many rows in the thing
        n_rows = math.ceil(len(self.exp_ids)/ncol)
        plt.subplots_adjust(bottom=0.18 + (n_rows-1)*0.06)
        es = []
        es.append(self.control_name)
        for e in self.exp_ids:
            es.append(self.experiment_names[e])
        fig.legend(pltHandles, es, ncol = ncol, loc='lower center', fontsize=10)


        # start the bit to make a nice fancy title.

        #make date strings to make date look prettier.
        sd = '{}'.format(self.control.startdate)
        ed = '{}'.format(self.control.enddate)
        startdate_fmt = sd[0:4]+'.'+sd[4:6]+'.'+sd[6:8]+' '+sd[8:10]+' UTC' 
        enddate_fmt = ed[0:4]+'.'+ed[4:6]+'.'+ed[6:8]+' '+ed[8:10]+' UTC'
        startdate,enddate = sd,ed

        n = self.cnt
        # Pick a title based on statistic
        if (statistic =='cor'):
            outtitle = "{} hPa Anomaly Correlation {} {}\n{} - {}, {} Forecasts".format(lev, var2Name[var], dom2Long[dom], startdate_fmt, enddate_fmt, n)
        elif(statistic =='rms'):
            outtitle = "{} hPa RMSE {} {}\n{} - {}, {} Forecasts".format(lev, var2Name[var], dom2Long[dom], startdate_fmt, enddate_fmt, n) 
        plt.suptitle(outtitle, fontsize=10)
        # done making fancy title.

        # output graphic.
        fn = '{}_z.{}.{}.{}.{}-{}.{}.verif-{}.{}'.format(statistic, var,lev, self.control.experiment, startdate, enddate, dom, verif1, verif2)
        print('Writing: ',fn + self.graphicOutput) 
        fig.savefig(fn + self.graphicOutput)

        #close fig (keep matplotlib from getting angry with too many open figures.)
        plt.close(fig)


    def make_significance_boxes(self, ax, xdata, ydata, xerror, yerror, facecolor='None',
                     edgecolor='None', alpha=0.5):
        """
        Basically taken from https://matplotlib.org/3.2.1/gallery/statistics/errorbars_and_boxes.html 
        Just changing defaults so it doesn't make things a shade of red and no real errorbar plot on top of it.. 
        """
        # Loop over data points; create box from errors at each point
        errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
                    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]

        # Create patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)

        # Add collection to axes
        ax.add_collection(pc)

        # put in "fake error bar" so it includes y-axis limits of the bars.
        artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          fmt='None', ecolor='None')
        return(artists)

if __name__ == "__main__":
    #a = plotter_2d()
    #a.plotter_2d()
    #a.plotter_2d(statistic='rms')

    a = plotter_1d()
    a.plotter_1d()
    a.plotter_1d(statistic='rms')
