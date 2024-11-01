def axvlines(ax=None, xs=[0, 1], ymin=0, ymax=1, **kwargs):
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    for x in xs:
        ax.axvline(x, ymin=ymin, ymax=ymax, **kwargs)

def mass_reconstruction_plot(matrix, events, total_reconst_mass, utotal_reconst_mass,unc,
                             event_columnname="Event_M",event_labels=['S', 'SN', 'SP','SL'],
                             savefig=True,imagepath='images/stacked_bar_daily_percentage_testM.png',
                             showplot=True):
    import matplotlib.pyplot as plt
    
    organic_mass_per = percentage_with_err(mass['organic_mass'], matrix['PM2.5'], uncertainty['uorganic_mass'], unc['PM2.5'])
    inorganic_ions_per = percentage_with_err(mass['inorganic_ions'], matrix['PM2.5'], uncertainty['uinorganic_ions'], unc['PM2.5'])
    geological_minerals_per = percentage_with_err(mass['geological_minerals'], matrix['PM2.5'], uncertainty['ugeological_minerals'], unc['PM2.5'])
    EC_per = percentage_with_err(mass['elemental_C'], matrix['PM2.5'], uncertainty['uelemental_C'], unc['PM2.5'])
    ssa_per = percentage_with_err( mass['salt'], matrix['PM2.5'], uncertainty['usalt'], unc['PM2.5'])
    others_per = ((mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3)/ total_reconst_mass * 100
    plt.style.use('seaborn-v0_8-paper')

    reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                              totalval=matrix['PM2.5'], utotalval=unc['PM2.5']) 
    fig, ax = plt.subplots(nrows=2, figsize=(7, 5), sharex=True, dpi=200)

    fig.suptitle('Mass reconstruction')
    # ax.set_title('Mass reconstructed')
    ax[0].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'],
                color='k', capsize=2, capthick=1, lw=1, marker='.', label='Gravimetric mass', zorder=1)
    ax[0].errorbar(matrix.index, total_reconst_mass, yerr=utotal_reconst_mass, color='red',
                capsize=2, capthick=1, lw=1, marker='.', label='Reconstructed mass', zorder=0)
    ax[0].set_ylabel('PM$_{2.5}$ (µg/m$^3$)')
    ax[0].plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels)) * 0, 'd',

            color='gray', label='Smoke events', zorder=3)
    # ax[0].plot(matrix.index, matrix['PM2.5'] - total_reconst_mass, '.-')
    # ax[0].plot(matrix.index, events[event_columnname].isin(event_labels), 'X')
    ax[0].legend()


    axvlines(ax=ax[0], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)
    axvlines(ax=ax[1], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)

    ax[1].bar(matrix.index.values, organic_mass_per['perc'].where(
        matrix['Na sol'].notna()).values, width,  label='OM')
    ax[1].bar(matrix.index.values, inorganic_ions_per['perc'].values,
            width,  bottom=organic_mass_per['perc'].values, label='II')
    ax[1].bar(matrix.index.values, geological_minerals_per['perc'].values, width,
            bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc']).values, label='GM')
    ax[1].bar(matrix.index.values, EC_per['perc'].values, width,
            bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']).values, label='EC')
    ax[1].bar(matrix.index.values, ssa_per['perc'].values, width,
            error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                        'ecolor': 'gray', 'marker': '.'},
            bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc'] + EC_per['perc']).values, label='SSA')
    ax[1].bar(matrix.index.values, others_per.values, width, yerr=reconst['uperc'],
            error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                        'ecolor': 'gray', 'marker': '.'},
            bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] +
                    geological_minerals_per['perc'] + EC_per['perc'] + ssa_per['perc']).values,
            label='Others')
    ax[1].axhline(100, linestyle=':', color='k')
    ax[1].axhline(100, linestyle=':', color='k')
    ax[1].axhspan(80, 120, alpha=0.2, color='y')
    ax[1].set_ylabel('Reconstructed mass (%)')
    ax[1].set_xlabel('Date')
    handles, labels = ax[1].get_legend_handles_labels()
    ax[1].legend(reversed(handles), reversed(labels), loc=1)
    fig.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    if savefig:
        fig.savefig(imagepath)
    if showplot:
        plt.show()