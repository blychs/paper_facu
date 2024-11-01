# Estima parametros lineales omoc 
for method in methods:
    resultNormal = linear_estimation_om_oc(matrix.where(
        events[event_columnname] == 'no'), method=method, ssa_as_Na=False, display_latex=True)
    omoc_noevent.append(resultNormal[0])
    resultEvent = linear_estimation_om_oc(matrix.where(
        events[event_columnname].isin(event_labels)), method=method,
        ssa_as_Na=False, display_latex=True)
    omoc_event.append(resultEvent[0])
    resultAll = linear_estimation_om_oc(matrix, method=method,
        ssa_as_Na=False, display_latex=True)
    omoc_all.append(resultAll[0])
    
    print(result.summary())

    # print(f"{method}")
    # print("No events")
    # print(resultNormal.summary())#.as_latex())
    # print("Events")
    # print(resultEvent.summary())#.as_latex())
    # print("All")
    # print(resultAll.summary())#.as_latex())
    
    print(method, " & No Events & ", f'{resultNormal[0]:.4g}', '&', f'{resultNormal[4]:.2g}', '&',
         f'{resultNormal[3]:.3g}', '&', f'{resultNormal[2]:.3g}', "& Events & ", f'{resultEvent[0]:.4g}', '&',
         f'{resultEvent[4]:.2g}', '&', f'{resultEvent[3]:.2g}', '&', 
         f'{resultEvent[2]:.3g}', " & All together & ", f'{resultAll[0]:.4g}', '&',
         f'{resultAll[4]:.2g}', '&',f'{resultAll[3]:.2g}', '&',
         f'{resultAll[2]:.3g}', '\\\\')

beta_omoc_noevent=np.round(np.mean(omoc_noevent),1)
beta_omoc_event=np.round(np.mean(omoc_event),1)
beta_omoc_all=np.round(np.mean(omoc_all),1)
print(beta_omoc_all,beta_omoc_event,beta_omoc_noevent)
