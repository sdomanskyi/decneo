import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from decneo.commonFunctions import getGenesOfPeak, movingAverageCentered, adjustTexts1D
from heatmap import Cluster,Avg3andUnique_Heatmap,MakeNoNa_RLavg3,makePath
import matplotlib.gridspec as gridspec
import seaborn as sns
from adjustText import adjust_text
from DBSCAN_fig import Group_DBSCAN, findGroupsDBSCAN

#lig_name = 'Fibroblast'
#rec_name = 'Endothelial'

def get_text_position(text, ax):
        x, y = text.get_position()
        return text.get_transform().transform((ax.convert_xunits(x), ax.convert_yunits(y)))
    
    


def adjustTexts1D_Y(texts, fig, ax, w = 'auto', direction = 'auto', maxIterations = 10**3, tolerance = 0.02):
    
    def get_text_position(text, ax):
        x, y = text.get_position()
        return text.get_transform().transform((ax.convert_xunits(x), ax.convert_yunits(y)))
    
    def set_text_position(text, x, y):
        return text.set_position(text.get_transform().inverted().transform((x, y)))

    extent = texts[0].get_window_extent(renderer=fig.canvas.get_renderer())
    
    if direction=='auto':
        direction = 'y' if extent.width > extent.height else 'x'
        
    if w == 'auto':      
        w = (extent.width if direction=='x' else extent.height) + 5
    
    orig_pos = [get_text_position(text, ax)[0 if direction=='x' else 1] for text in texts]
    
    for i, text in enumerate(texts):
        x, y = get_text_position(text, ax)
        set_text_position(text, x + i*w*(1 if direction=='x' else 0), y + i*w*(0. if direction=='x' else 1))
    
    objs = []
    for i_iter in range(maxIterations):
        curr_pos = [get_text_position(text, ax)[0 if direction=='x' else 1] for text in texts]
        
        obj = 0
        for i, (opos, text) in enumerate(zip(orig_pos, texts)):
            x, y = get_text_position(text, ax)
            p, q = (y, x) if direction=='x' else (x, y)
            
            cx = curr_pos[np.argsort(np.abs(curr_pos - q))[1]]
            ov, ovdel = max(0, w + min(q, cx) - max(q, cx)), max(0, w + min(q, cx+0.001) - max(q, cx+0.001))
            
            delta1 = 0.1*ov*np.sign(ovdel - ov)
            delta2 = np.sign(q - opos) * np.abs(q - opos)/25.

            delta = delta1 - delta2 * ((w-ov)/w)**2.5
            
            set_text_position(text, *((q + delta, p) if direction=='x' else ((-2*p) , q + (delta))))
            
            obj += 5.*ov + np.abs(q - opos)
        objs.append(obj)
        
        try:
            if abs(np.mean(objs[-20:]) - np.mean(objs[-40:-20])) < tolerance: 
                break
        except:
            pass
        
    return



def Make_Combined_Graph(lig_name,rec_name, DBSCAN = False,Comberon_Box = False):
    fname = 'intermediate_choroid/'+lig_name+'_'+rec_name+'_RLavg3mean_scaled_nona.xlsx'

    cutoff = 0.0
    lig = pd.read_excel('results 09 17 2021 '+str(cutoff)+'/ligands/'+lig_name+'. dendrogram-heatmap-correlation-data.xlsx', index_col = 0)
    rec = pd.read_excel('results 09 17 2021 '+str(cutoff)+'/receptors/'+rec_name+'. dendrogram-heatmap-correlation-data.xlsx', index_col = 0)

    fig = plt.figure()
    #ax1 = plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=2) # topleft


    
    # 0.02
    width_ratios = [0.04,0.02,0.04,1] #y,x,heatmap
    height_ratios = [0.1,0.1,0.2,1,1,1,1,1,1,1]


    fig = plt.figure(figsize = (15,16),tight_layout = True)
    fig.suptitle(lig_name+' x '+rec_name +' RLavg3 combo scaled score', size = 20)


    gs = fig.add_gridspec(10,4, width_ratios = width_ratios, height_ratios = height_ratios)
    gs.update(wspace=0.025,hspace = 0.05)
    
    ax1 = fig.add_subplot(gs[1, 3:4]) #receptors cluster
    ax2 = fig.add_subplot(gs[2, 3:4],sharex = ax1) #receptors
    ax3 = fig.add_subplot(gs[3:10,1]) #ligands cluster
    ax4 = fig.add_subplot(gs[3:10,2],sharey = ax3) #ligands

    main_ax = fig.add_subplot(gs[3:10,3:4], sharex = ax1, sharey = ax3) #heatmap

    #ax3.set_aspect(aspect = 1)
    #ax4.set_aspect(aspect = 1)
    
    width = main_ax.get_xticks()[1] - main_ax.get_xticks()[0]
    
    
    lig_combo3avgs_smooth = lig['Avg Combination of measures']
    rec_combo3avgs_smooth = rec['Avg Combination of measures']

    lig_combo3avgs = lig['Combination of measures']
    rec_combo3avgs = rec['Combination of measures']
    
    ax2.set_ylim(0,round(max(rec_combo3avgs),3))
    ax4.set_xlim(0,round(max(lig_combo3avgs),3))


    #cbar_ax = plt.subplot2grid((14,14), (1,13), rowspan = 13)



    fig.tight_layout()




    #ax3.plot(z['Avg Combination of measures'])


    lig_ncluster = len(lig['cluster'].unique().tolist())
    rec_ncluster = len(rec['cluster'].unique().tolist())
    lig_peaks = getGenesOfPeak(pd.Series(lig_combo3avgs/lig_combo3avgs.max()))
    rec_peaks = getGenesOfPeak(pd.Series(rec_combo3avgs/rec_combo3avgs.max()))
    
    

    lig_bw = ax4.get_ylim()[0]/lig_ncluster
    rec_bw = ax3.get_xlim()[1]/rec_ncluster






    #data_avg = movingAverageCentered(np.nan_to_num(data), 10)
    lig_vals = lig['Combination of measures'].values
    y_vals = np.array([b for b in range(0,len(lig_vals))])

    lig_vals_smooth = lig['Avg Combination of measures'].values



    main_ax.set_xticks([])
    main_ax.set_yticks([])
    main_ax.axes.xaxis.set_visible(False)
    main_ax.axes.yaxis.set_visible(False)
    
    #ax1.axes.xaxis.set_visible(False)
    ax1.axes.yaxis.set_visible(False)
    ax1.xaxis.tick_top()
    ax1.axes.xaxis.set_ticks([])
    ax1.axes.xaxis.set_ticklabels([])
    
    ax2.axes.xaxis.set_visible(False)
    ax2.axes.yaxis.set_visible(False)
    
    ax3.axes.xaxis.set_visible(False)
    #ax3.axes.yaxis.set_visible(False)
    
    ax4.axes.xaxis.set_visible(False)
    ax4.axes.yaxis.set_visible(False)
    if DBSCAN:
        z = Group_DBSCAN(fname,fig = fig,ax = main_ax,eps = 10, min_samples = 50, plot = True, findMaxVal = False, getTop = 1)
    else:
        z = Cluster(fname, ax = main_ax,cbar_ax = None, fig = fig, plot = True, minDist = 30, minSizeClust = 50, minSizeBounds = 50, getTop = 1, findMax = False)
    
    

    
    ax2.plot(rec_combo3avgs,zorder = 1)
    ax2.plot(rec_combo3avgs_smooth,zorder = 4)

    ax2.fill_between(rec_combo3avgs.index,rec_combo3avgs.values,color = 'darkblue', zorder = 3)

    ax4.plot(lig_vals,y_vals,zorder = 1)
    ax4.plot(lig_vals_smooth,y_vals,zorder = 4)
    ax4.fill_betweenx(lig_combo3avgs.index,lig_vals,color = 'darkblue', zorder = 3)

    peak_recs_index = [rec.index.get_loc(z) for z in rec_peaks]
    peak_ligs_index = [lig.index.get_loc(z) for z in lig_peaks]

    max_lig, min_lig = max(peak_ligs_index), min(peak_ligs_index)
    max_rec, min_rec = max(peak_recs_index), min(peak_recs_index)

    ax2.axvspan(min_rec, max_rec, facecolor='crimson',zorder = 2)
    ax4.axhspan(min_lig, max_lig, facecolor='crimson',zorder = 2)
    
    if Comberon_Box:
        main_ax.plot([min_rec,min_rec], [min_lig, max_lig], marker = 'o',ls = '--',c = 'crimson')
        main_ax.plot([max_rec,max_rec], [min_lig,max_lig], marker = 'o',ls = '--',c = 'crimson')
        main_ax.plot([min_rec,max_rec],[min_lig, min_lig], marker = 'o',ls = '--',c = 'crimson')
        main_ax.plot([min_rec, max_rec], [max_lig,max_lig], marker = 'o',ls = '--',c = 'crimson')
        Comberon_dict = {}
        #print(min_rec)
        #print(max_lig)
        comb_ligs = [lig.index[z] for z in range(min_lig,max_lig)]
        comb_recs = [rec.index[z] for z in range(min_rec,max_rec)]
        if len(comb_ligs) > len(comb_recs):
           while len(comb_ligs) > len(comb_recs):
               comb_recs.append(np.nan)
        else:
            while len(comb_recs) > len(comb_ligs):
                comb_ligs.append(np.nan)
        
        Comberon_dict['ligs'] = comb_ligs
        Comberon_dict['recs'] = comb_recs
        df = pd.DataFrame.from_dict(Comberon_dict)
        df.to_excel('intermediate_choroid/'+lig_name+'_'+rec_name+'_comberon.xlsx', index = False)
    
    
    lig_cluster = lig[['cluster']]
    
    rec_cluster = rec[['cluster']]
    
    
    
    
    unique_cluster_lig = np.unique(lig_cluster['cluster']).tolist()
    unique_cluster_rec = np.unique(rec_cluster['cluster']).tolist()
    top_n_cluster = 1
    if len(unique_cluster_lig) > len(unique_cluster_rec):
        top_n_cluster = len(unique_cluster_lig)
    else:
        top_n_cluster = len(unique_cluster_rec)
        
    possible_colors = sns.color_palette('Paired', top_n_cluster)
    
    for cluster in unique_cluster_rec:
        rec_subset = rec_cluster[rec_cluster['cluster'] == cluster]
        ax1.bar(rec_subset.index.values, height = 1, width = 0.02,edgecolor = possible_colors[cluster])
        ax1.axes.xaxis.set_ticks([])
        ax1.axes.xaxis.set_ticklabels([])
    for cluster in unique_cluster_lig:
        lig_subset = lig_cluster[lig_cluster['cluster'] == cluster]
        ax3.barh(lig_subset.index.values, width = 0.02, height = 1, align = 'center', edgecolor = possible_colors[cluster])
        
        ax3.axes.yaxis.set_ticks([])
        ax3.axes.yaxis.set_ticklabels([])
    
    if DBSCAN:
        ax3.invert_yaxis() #turn this on if DBSCAN method otherwise off
        
    
    interesting_choroid_lig = ['ADIPOQ', 'ANGPT1', 'ANGPT2', 'ANGPTL3', 'BMP2', 'BMP7', 'C3', 'C4B', 'DLL1',
                              'FN1', 'IL13', 'IL4', 'JAG1', 'JAG2', 'LIF', 'OSM', 'S1P', 'SEMA3A', 'SEMA3C',
                              'SEMA3E', 'SEMA4A', 'SLIT2', 'TGFB1', 'TGFB2', 'TGFB3', 'VEGFA'] # 26
    
    interesting_choroid_rec = ['KDR','FLT1','FLT4','NRP1','NRP2','FGFR1','FGFR2','FGFR3','CXCR2','ROBO1',
             'ROBO4','ENG','PDGFRA','PDGFRB','TEK','KIT','MET','CLEC14A', # stimulators
             'CD36','CD47','VLDLR','PLXND1']
    
    lig_interesting_choroid_in = lig.loc[lig.index.intersection(interesting_choroid_lig),:].index.values.tolist()
    lig_indicies = []
    for z in lig_interesting_choroid_in:
        lig_indicies.append(lig.index.get_loc(z))
    rec_interesting_choroid_in = rec.loc[rec.index.intersection(interesting_choroid_rec),:].index.values.tolist()
    rec_indicies = []
    for z in rec_interesting_choroid_in:
        rec_indicies.append(rec.index.get_loc(z))
        
    
    
    rec_texts = [ax1.text(rec_indicies[ngene],4.5,rec_interesting_choroid_in[ngene],fontsize = 12, rotation = 90) for ngene in range(len(rec_interesting_choroid_in))]
    lig_texts = [ax3.text(-0.25,lig_indicies[ngene],lig_interesting_choroid_in[ngene], fontsize = 12) for ngene in range(len(lig_interesting_choroid_in))]
    
    rec_an = adjustTexts1D(rec_texts,fig,ax=ax1)
    if lig_name in ['Pericyte','Endothelial']:
        lig_an = adjust_text(lig_texts, ax = ax3, va = 'center', ha = 'right', autoalign = 'y', lim = 400, force_text = (3,4))
    else:
        lig_an = adjust_text(lig_texts, ax = ax3, va = 'center', ha = 'right', autoalign = 'y', lim = 400, expand_text = (0.5,0.9))
    v_rec = float(ax1.get_ylim()[1])
    orig_pos_rec = [rec_indicies[ngene] for ngene in range(len(rec_interesting_choroid_in))]
    
    h_lig = float(ax3.get_xlim()[1])
    orig_pos_lig = [lig_indicies[ngene] for ngene in range(len(lig_interesting_choroid_in))]
    

    
    for text, opos in zip(rec_texts,orig_pos_rec):
        text._y = 3 * v_rec
        ax1.plot([text._x, opos], [text._y,1], color = text._color, lw = 0.5, clip_on = False)
        
        
    #lig_texts = [ele for ele in reversed(lig_texts)]
    for text, opos in zip(lig_texts,orig_pos_lig):
        text._x = -2.5 * h_lig
        ax3.plot([0.,text._x ], [opos,text._y], color = text._color, lw = 0.5, clip_on = False)
    
    ax1.set_facecolor('white')
    ax3.set_facecolor('white')
    
    """
    for ngene in range(len(rec_interesting_choroid_in)):
        an_rec = ax1.annotate(text = rec_interesting_choroid_in[ngene],xy = (rec_indicies[ngene],1),xytext = (rec_indicies[ngene], 2.5), xycoords = 'data',textcoords = 'data',
                     arrowprops = dict(arrowstyle = '-', color = 'green'), fontsize = 20,zorder = 4)
        
    for ngene in range(len(lig_interesting_choroid_in)):
        an_lig = ax3.annotate(text = lig_interesting_choroid_in[ngene],xy = (0,lig_indicies[ngene]),xytext = (-0.2,lig_indicies[ngene]), xycoords = 'data',textcoords='data',
                     arrowprops = dict(arrowstyle = '-', color = 'green'), fontsize = 20,zorder = 4)
    """
        
    
    
        
    
    

    #ax3.plot(rec['Combination of measures'])
    

    if DBSCAN:
        fpath = 'out_choroid/heatmaps/DBSCAN/RL_avg3/'
        makePath(fpath)
        fig.savefig(fpath+lig_name+'_'+rec_name+'_RLavg3_scaled.png',dpi = 1000)
    else:
        fpath = 'out_choroid/heatmaps/RL_avg3/'
        makePath(fpath)
        fig.savefig(fpath +lig_name+'_'+rec_name+'_RLavg3_scaled.png',dpi = 1000)
    
    plt.clf()
    
    return True


lig_names = ['Fibroblast'] # 'Fibroblast','Macrophage','SMC','Pericyte','Endothelial'
rec_names = ['Endothelial']
thresh = 0.5
for lig_name in lig_names:
    for rec_name in rec_names:
        print(lig_name + ',' + rec_name)
        start = Avg3andUnique_Heatmap(ligand = lig_name,receptor = rec_name, mode = 'sum', plot = False)
        start1 = MakeNoNa_RLavg3('intermediate_choroid/'+lig_name+'_'+rec_name+'_RLavg3mean_scaled.xlsx',thresh = thresh)
        b = Make_Combined_Graph(lig_name = lig_name, rec_name = rec_name, DBSCAN = True)
        a = Make_Combined_Graph(lig_name = lig_name, rec_name = rec_name, DBSCAN = False,Comberon_Box = False)
        

