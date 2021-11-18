import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.ticker import PercentFormatter
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import cmocean
import csv
import copy
import os
from pathlib import Path
from scipy import stats
import statsmodels.api as sm 

import config as plot_config
#################################################

# * * *
#################################################
def make_con_csv(Surge_ID, sigma, mu, kappa, gamma, day, Data, config):
    ####
    if config.Conv_controls['Valid_Rays']== True:
        csv_name= config.g_TORUS_directory+'VRays_Convergence_stats.csv'
    else:
        csv_name= config.g_TORUS_directory+'Convergence_stats.csv'
    Does_csv_exist=os.path.isfile(csv_name)
    if Does_csv_exist == False: 
        # making csv file if it does not exist 
        with open(csv_name, 'w') as csvfile: 
            # field names 
            fields = ['Day','RName', 'Scan_time', 'Surge_ID', 'Sigma', 'Mu', 'Kappa','Gamma'] 
            
            # creating a csv writer object 
            csvwriter = csv.writer(csvfile) 
            # writing the fields 
            csvwriter.writerow(fields) 
    ###
    # Now append the new row of data points
    with open(csv_name, 'a') as csvfile: 
        # writing the data rows 
        csvwriter = csv.writer(csvfile) 
        Rname= Data['P_Radar'].name
        rows= [[day, Rname,Data['P_Radar'].Scan_time, Surge_ID, sigma, mu, kappa, gamma]]
        csvwriter.writerows(rows)

#################################################
def csvstats_scatterplts(config, compday1, compday2):
    print('made it to here')
    if config.Conv_controls['Valid_Rays']== True:
        csv_name= config.g_TORUS_directory+'VRays_Convergence_stats.csv'
    else:
        csv_name= config.g_TORUS_directory+'Convergence_stats.csv'
    df = pd.read_csv(csv_name)

    #  betterday=df[df['Day']== 20190524]
    #  subdf=df.loc[(df['Day']==20190524).all() or (df['Day']== 20190617).all()]
    #  subdf=df[df['Day']==20190524 or df['Day']== 20190617]
    #  compday1=20190524
    #  compday2=20190615
    subdf1=df[df['Day']==compday1].reset_index()# or df['Day']== 20190617]
    subdf2=df[df['Day']==compday2].reset_index()# or df['Day']== 20190617]

    kappa_outliers1=subdf1[subdf1['Kappa']>100]
    gamma_outliers1=subdf1[subdf1['Gamma']>20]
    kappa_outliers2=subdf2[subdf2['Kappa']>100]
    gamma_outliers2=subdf2[subdf2['Gamma']>20]
    num_of_Koutliers=len(kappa_outliers1)+len(kappa_outliers2)
    num_of_Goutliers=len(gamma_outliers1)+len(gamma_outliers2)
    #  print(subdf1)
    #  print('ooooooooooooooooo')
    #  print(subdf2)
    fig = plt.figure(figsize=(9,10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=.1)


    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    stats=['Sigma', 'Mu', 'Kappa', 'Gamma']
    axis=[ax0, ax1, ax2, ax3]

    for ax, stat in zip(axis, stats):
        ax.scatter(subdf1.index, subdf1.loc[:,stat], color='g', label=compday1)
        ax.scatter(subdf2.index, subdf2.loc[:,stat], color='r', label=compday2)

        m1, b1= np.polyfit(subdf1.index, subdf1.loc[:, stat], 1)
        m2, b2= np.polyfit(subdf2.index, subdf2.loc[:, stat], 1)
        ax.plot(subdf1.index, m1*subdf1.index+b1, color='g', zorder=7, alpha=.8)
        ax.plot(subdf2.index, m2*subdf2.index+b2, color='r', zorder=7, alpha=.8)

    #  ax0.scatter(df.index, df.loc[:,'Sigma'])
    #  ax1.scatter(df.index, df.loc[:,'Mu'])
    #  ax2.scatter(df.index, df.loc[:,'Kappa'])
    #  ax3.scatter(df.index, df.loc[:,'Gamma'])
    #  ax0.scatter(subdf1.index, subdf1.loc[:,'Sigma'], color='g', label=compday1)
    #  ax1.scatter(subdf1.index, subdf1.loc[:,'Mu'], color='g')
    #  ax2.scatter(subdf1.index, subdf1.loc[:,'Kappa'], color='g')
    #  ax3.scatter(subdf1.index, subdf1.loc[:,'Gamma'], color='g')
    #  ax0.scatter(subdf2.index, subdf2.loc[:,'Sigma'], color='r', label=compday2)
    #  ax1.scatter(subdf2.index, subdf2.loc[:,'Mu'], color='r')
    #  ax2.scatter(subdf2.index, subdf2.loc[:,'Kappa'], color='r')
    #  ax3.scatter(subdf2.index, subdf2.loc[:,'Gamma'], color='r')
#
    #  m1, b1= np.polyfit(subdf1.index, subdf1.loc[:, 'Sigma'], 1)
    #  ax0.plot(subdf1.index, m1*subdf1.index+b1)
#
    ax0.set_ylabel('Std Dev ('+r'$\sigma$'+')')
    ax1.set_ylabel('Mean ('+r'$\mu$'+')')
    ax2.set_ylabel('Kurtosis ('+r'$\kappa$'+')')
    ax3.set_ylabel('Skew ('+r'$\gamma$'+')')
    ax0.legend()
    ax2.set_ylim(top=100)
    ax3.set_ylim(top=20)
    ax2.text(len(subdf1),95,str(num_of_Koutliers))
    ax3.text(len(subdf1),15,str(num_of_Goutliers))

    if config.Conv_controls['Valid_Rays']== True:
        plt.suptitle(str(compday1)+' vs '+str(compday2)+'VRays',ha='center')
        outdir = config.g_TORUS_directory+'VRays_Convergence_stats_plt_'+str(compday1)+'vs'+str(compday2)+'.png'
    else:
        plt.suptitle(str(compday1)+' vs '+str(compday2),ha='center')
        outdir = config.g_TORUS_directory+'Convergence_stats_plt_'+str(compday1)+'vs'+str(compday2)+'.png'
    #  plt.tight_layout()
    plt.savefig(outdir)
    plt.close()

#################################################
def Conv_1dplt(Z, X, Y, Surge_ID, day, Data, config, tilt_ang, valid_rays, Surge_stats_dict):
    print('Made it to Conv_1dplt')
    cmaps = ['PuOr', 'PiYG', cmocean.cm.balance, cmocean.cm.balance,cmocean.cm.dense]
    N=[colors.Normalize(-.03, .03), colors.Normalize(-20, 20), colors.Normalize(-30, 30), colors.Normalize(-30, 30), None]
    lab=['Conv', 'Vel Diff', 'Minus Rbins Vel', 'Plus Rbins Vel', 'Dist']
    under_c=['#FF3131','magenta','aqua', 'aqua', 'aqua']
    over_c=['magenta','lime','magenta', 'magenta', 'magenta']
    
    #  fig, axs = plt.subplots(3, sharex=True)
    fig = plt.figure(figsize=(9,10))
    gs = gridspec.GridSpec(5, 1, figure=fig, height_ratios=[3,3,3,3,1], hspace=.1)

    for row in range(5):
        if config.Conv_controls['Valid_Rays']== True:
            for ray in range(np.shape(Z[row])[0]):
                if valid_rays[ray]==False:
                    Z[row][ray,:]=np.nan
        ax = fig.add_subplot(gs[row, :])

        colormap = copy.copy(matplotlib.cm.get_cmap(name=cmaps[row]))
        colormap.set_under(under_c[row])
        colormap.set_over(over_c[row])

        #  print('X: {}, Y: {}, Z: {}, Zrow: {}, N: {}, Nrow: {}, CM: {}'.format(np.shape(X),np.shape(Y),np.shape(Z),np.shape(Z[row]),np.shape(N),np.shape(N[row]),colormap))
        pcm = ax.pcolormesh(X, Y, Z[row], norm=N[row], cmap=colormap)
        ax.set_ylabel(lab[row])
        ax.yaxis.set_label_position("right")

        if row != 4:
            ax.tick_params(labelbottom=False)

        fig.colorbar(pcm, ax=ax, extend='both')

    if config.Conv_controls['Valid_Rays']== True:
        sweep = 'Vrays_'+str(Data['P_Radar'].swp)
    else:
        sweep = str(Data['P_Radar'].swp)
    plt.suptitle(Data['P_Radar'].site_name+' '+str(sweep)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str+' '+str(Surge_ID)+' ('+str(tilt_ang[0])+')', ha='center')

    #path to file 
    plt_dir=Data['P_Radar'].dir_name+'/clicker/convergence/'+str(sweep)+'_deg/'
    #add the file name
    output_name = Data['P_Radar'].site_name+'_'+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+str(Surge_ID)+'.png'
    #This is the directory path for the output file
    outdir = config.g_TORUS_directory+day+'/plots/'+plt_dir
    if not os.path.exists(outdir): Path(outdir).mkdir(parents=True, exist_ok=True)
    output_path_plus_name = outdir + output_name
    print(output_path_plus_name)
    plt.tight_layout()
    plt.savefig(output_path_plus_name)
    plt.close()


#################################################
def Conv_histograms(Conv, Surge_ID, day, Data, config, tilt_ang, valid_rays, Surge_stats_dict):
    print('Made it to Conv_histograms')
    ## set up the data to put into histogram
    if config.Conv_controls['Valid_Rays']== True:
        for ray in range(np.shape(Conv)[0]):
            if valid_rays[ray]==False:
                Conv[ray,:]=np.nan
    Conv=np.ndarray.flatten(Conv)
    Conv_nonans=Conv[~np.isnan(Conv)]
    if len(Conv_nonans)==0:
        return 
    sigma=np.std(Conv_nonans)
    mu=np.mean(Conv_nonans)
    kappa=stats.kurtosis(Conv_nonans)
    gamma=stats.skew(Conv_nonans)

    if config.Conv_controls['add_to_con_csv']==True:
        make_con_csv(Surge_ID, sigma, mu, kappa, gamma, day, Data, config)
    #  csvstats_scatterplts(config)


    ## set up the histogram bins
    cmin=-.03 #np.nanmin(Conv)
    cmax=.03 #np.nanmax(Conv)
    bwidths=[.001, .0025, .005, .01]
    #  interval=.0025
    
    fig = plt.figure(figsize=(15,6))
    gs = gridspec.GridSpec(1, 5, figure=fig, width_ratios=[3,3,3,3,1])
    ###

    ax = fig.add_subplot(gs[:, -1])

    ax.text(0, .98, Data['P_Radar'].site_name+' '+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+' '+str(Surge_ID), ha="left", va="top", weight='bold', size=15)
    ax.text(0, .88, 'Sigma: '+str("%.5f" %sigma), ha="left", va="top")
    ax.text(0, .78, 'Mu: '+str("%.5f" %mu), ha="left", va="top")
    ax.text(0, .68, 'Kappa: '+str("%.5f" %kappa), ha="left", va="top")
    ax.text(0, .58, 'Gamma: '+str("%.5f" %gamma), ha="left", va="top")
    #  ax.text(0, .48, 'tilt_ang: '+str(tilt_ang), ha="left", va="top")
    ax._frameon=False
    ax.tick_params(which='both', axis='both', left=False, bottom=False, labelleft=False, labelbottom=False)





    ###
    for pltnum, interval in enumerate(bwidths):
        ax = fig.add_subplot(gs[:, pltnum])

        #the list of bins including an extra on each end to capture any outliers (plus the fact that arange is not inclusive)
        bins_list=np.arange(cmin-interval, cmax+(2*interval), interval)
        norm=plt.Normalize(vmin=cmin,vmax=cmax)

        ## make histogram
        #  fig,ax = plt.subplots(figsize=(8,9))
        n, bins, patches=plt.hist(np.clip(Conv, bins_list[0], bins_list[-1]), edgecolor='black', bins=bins_list, zorder=10)
        #  plt.hist(Conv, bins=bins_list, density=True)
        
        ## apply coloring that matches the 2d conver fig
        colormap = copy.copy(matplotlib.cm.get_cmap(name='PuOr'))
        c=colormap(norm(bins_list))
        for i in range(len(patches)):
            if i == len(patches)-1: color='magenta'
            elif i == 0: color='#FF3131'
            else: color=c[i]
            patches[i].set_facecolor(color)

        ## set up axis details
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=len(Conv_nonans)))
        #  ax.set_xlim(left=-.03, right=.03)
        #find where 30% is for the total number of radarbins
        maxper_rbins=(40 * len(Conv_nonans))/100
        incper_rbins=(5 * len(Conv_nonans))/100
        minor_incper_rbins= incper_rbins/2
        maj_ticks=np.arange(0, maxper_rbins+incper_rbins, incper_rbins) 

        ax.yaxis.set_major_locator(ticker.FixedLocator(maj_ticks))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.set_ylim(top=maxper_rbins)
        ax.grid(which='major', axis='y', zorder=7, alpha=.8)
        ax.grid(which='minor', axis='y', zorder=7, alpha=.8, linestyle=':')
        ax.axvline(mu, linestyle=':', alpha=.8, zorder=10)

        if pltnum == 0:
            ax.set_ylabel('Percentage of Radarbins')
        ax.text(-.03, 1, 'bwidth='+str(interval))
        #  ax.text(1, 0, 'bwidth='+str(interval))
        #  ax.text(0, 0, 'bwidth='+str(interval), ha="left", va="top")
        ax.label_outer()

    #  xlabels=bins_list[1:].astype(str)
    #  xlabels[-1] += '+'
    #  ax.set_xticklabels(xlabels)

    #  plt.title(Data['P_Radar'].site_name+' '+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+' '+Surge_ID)
            #  + ': bwidth ='+str(interval)+'\n'+r'$\mu=$'+str('%.3E'%mu)+', '+r'$\sigma=$'+str('%.3E'%sigma)+',\n '+r'$\kappa=$'+str('%.3E'%kappa)+', '+r'$\gamma=$'+str('%.3E'%gamma), size=10)
    ax.set_xlabel('Convergence /s')

    plt_dir=Data['P_Radar'].dir_name+'/clicker/convergence/'+str(Data['P_Radar'].swp)+'_deg/'
    outdir = config.g_TORUS_directory+day+'/plots/'+plt_dir+'histogram/'
    if not os.path.exists(outdir): Path(outdir).mkdir(parents=True, exist_ok=True)
    output_name = 'hist_'+Data['P_Radar'].site_name+'_'+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+str(Surge_ID)+'.png'
    output_path_plus_name = outdir + output_name
    print(output_path_plus_name)
    plt.tight_layout()
    plt.savefig(output_path_plus_name)
    plt.close()


    #  QQplt(Conv_nonans, Surge_ID, day, Data, config)


def QQplt(obs, Surge_ID, day, Data, config):
    #standardize ob
    z=(obs-np.mean(obs))/np.std(obs)

    csv_name = config.g_TORUS_directory+'/z_default.csv'
    if Data['P_Radar'].Scan_time.strftime('%m%d_%H%M') =='0524_2033':
        Does_csv_exist = os.path.isfile(csv_name)
        if Does_csv_exist == False: 
            # making csv file if it does not exist 
            with open(csv_name, 'w') as csvfile: 
                # field names 
                fields = ['Day','Radar','Time','Surge_ID','z_default'] 
                
                # creating a csv writer object 
                csvwriter = csv.writer(csvfile) 
                # writing the fields 
                csvwriter.writerow(fields) 
        ###
        # Now append the new row of data points
        with open(csv_name, 'a') as csvfile: 
            # writing the data rows 
            csvwriter = csv.writer(csvfile) 
            rows = [[day, Data['P_Radar'].site_name, Data['P_Radar'].Scan_time.strftime('%H%M'), Surge_ID, z]]
            csvwriter.writerows(rows)

    df = pd.read_csv(csv_name)
    zdefault = df.loc[:,'z_default'].to_numpy()
    print(zdefault)
    print(type(zdefault))
    
    fig=sm.qqplot(obs, line='45')

    #  stats.probplot(z, dist='norm', plot=plt)
    #  statsmodels.graphics.gofplots.ProbPlot.qqplot()
    #  countabove= len([i for i in z if i > 10])
    #  countbelow= len([i for i in z if i < 10])
    #  count=countabove+countbelow
    #  plt.ylim(bottom=-10, top=10)
    #  plt.title('Nomal QQ plot\n'+Data['P_Radar'].site_name+' '+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+' '+str(Surge_ID)+' ('+str(count)+')')
    plt.title('Nomal QQ plot\n'+Data['P_Radar'].site_name+' '+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+' '+str(Surge_ID))
    
    plt_dir=Data['P_Radar'].dir_name+'/clicker/convergence/'+str(Data['P_Radar'].swp)+'_deg/'
    outdir = config.g_TORUS_directory+day+'/plots/'+plt_dir+'QQ/'
    if not os.path.exists(outdir): Path(outdir).mkdir(parents=True, exist_ok=True)
    output_name = 'QQ_'+Data['P_Radar'].site_name+'_'+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+str(Surge_ID)+'.png'
    output_path_plus_name = outdir + output_name
    print(output_path_plus_name)
    plt.tight_layout()
    plt.savefig(output_path_plus_name)
    plt.close()


######################################################
if plot_config.Conv_controls['make_scatterplts']==True:
    for day2 in plot_config.day_list:
        csvstats_scatterplts(plot_config, compday1=20190524, compday2=int(day2))

