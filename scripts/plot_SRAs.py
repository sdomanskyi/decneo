def plotDf(df, figName):
    fig = plt.figure(figsize=(8,10))

    n_clusters = 10
    ax = fig.add_axes([0.35, 0.7, 0.55, 0.15], frame_on=False)
    Z = hierarchy.linkage(np.nan_to_num(df.values, nan=df.values.max()), method='ward', optimal_ordering=True)
    origLineWidth = matplotlib.rcParams['lines.linewidth']
    matplotlib.rcParams['lines.linewidth'] = 0.5
    D = hierarchy.dendrogram(Z, ax=ax, color_threshold=0, above_threshold_color='k', orientation='top')
    matplotlib.rcParams['lines.linewidth'] = origLineWidth
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])

    plt.title(figName)

    tissues = [c.split(' Homo sapiens ')[0] if 'Homo sapiens' in c else c.split(' Mus musculus ')[0] for c in df.columns]
    dcolors = {c[1]:c[0] for c in list(enumerate(np.unique(tissues)))}
    colors = {k: cm.jet(v/len(dcolors)) for k, v in dcolors.items()}
    colors = [colors[t] for t in tissues]
    colors = np.array(colors)[D['leaves']]

    ax = fig.add_axes([0.35, 0.3, 0.55, 0.4], frame_on=False)
    cmap = plt.cm.hot
    cmap.set_bad('grey')
    sdata = np.ma.array(df.values[:,D['leaves']][D['leaves'],:], mask=np.isnan(df.values))
    im = ax.imshow(sdata, cmap=cmap, aspect='auto', interpolation='None', extent=(-0.5, df.shape[0] - 0.5, df.shape[1] - 0.5, -0.5))
    ax.set_xticks(range(len(df.columns)))
    ax.set_yticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns.values[D['leaves']])
    ax.set_yticklabels(df.columns.values[D['leaves']])
    ax.tick_params(axis='x', labelsize=8, width=0.25, length=1, rotation=90)
    ax.tick_params(axis='y', labelsize=8, width=0.25, length=1, rotation=0)

    for i, tick in enumerate(plt.gca().get_xticklabels()):
        tick.set_color(colors[i])

    for i, tick in enumerate(plt.gca().get_yticklabels()):
        tick.set_color(colors[i])

    ax = fig.add_axes([0.9, 0.5, 0.025, 0.25], frame_on=False)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    clb = fig.colorbar(im, ax=ax, fraction=0.4, label='Bootstrap counts corr.')
    clb.ax.tick_params(labelsize=4)

    plt.savefig(figName, dpi=300)
    plt.clf()

    return