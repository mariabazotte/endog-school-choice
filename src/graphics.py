import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

# Define file and save paths
file_path = '../results/compiled/results_paper_final.xlsx'
file_path_small = '../results/compiled/results_paper_small.xlsx'
save_path = '../results/compiled/graphics/'

#############################################################################################################
#                          Deterministic (average scenario) vs. Stochastic problems                         #
#############################################################################################################

# Information
sheet_names = ['vss_tr0', 'vss_tr1', 'vss_tr3']
label_names = {'vss_tr0':'UM','vss_tr1':'CEUM','vss_tr3':'IEUM'}
# colors = {'UM':'tab:blue','CEUM':'tab:green','IEUM':'tab:orange'}
colors = {'UM': '#0072B2','CEUM': '#009E73','IEUM': '#E69F00'}
style  = {'UM':'-','CEUM':'--','IEUM':'-.'}
hatch_style = {'UM':'/','CEUM':'x','IEUM':'o'}
det_stoch_path = save_path + 'det_stoch/'

def plot_metric_by_budget(y_column, y_column_name):
    """
    Plots the average and quantiles of a metric (y_column) by budget for multiple types of student behavior, and saves the figure.
    These metrics correspond to information comparing the stochastic vs. deterministic version of the problem.
    
    Parameters:
    - y_column: Name of the column to use for the y-axis. ('vss','nb_enter','nb_improv','nb_total')
    - y_column_name: Name to display on the y-axis 
    """
    
    all_maxranks = []
    plt.figure(figsize=(16, 10))
    for sheet in sheet_names:
        df = pd.read_excel(file_path, sheet_name=sheet)
        if y_column not in df.columns:
            print(f"Column '{y_column}' not found in sheet '{sheet}' — skipping.")
            continue
        grouped = df.groupby('maxseats', as_index=False)[y_column].mean()
        grouped = grouped.sort_values(by='maxseats')
        all_maxseats = sorted(df['maxseats'].unique())
        
        name = label_names.get(sheet, None)
        plt.plot(grouped['maxseats'], grouped[y_column], label = rf'$\texttt{{{name}}}$', color=colors[name], linestyle=style[name], marker='o', linewidth=5)

    plt.xlabel(r'Budget $B$')
    plt.ylabel(f'{y_column_name}')
    plt.ylim(bottom=0) 
    plt.xticks(all_maxseats) 
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(det_stoch_path + f"{y_column}.png", dpi=300)
    plt.close()
    
    #########################################
    
    plt.figure(figsize=(16, 10))
    
    aux = {r'$\texttt{UM}$':'UM',r'$\texttt{CEUM}$':'CEUM',r'$\texttt{IEUM}$':'IEUM'}
    aux_hatch = ['UM','CEUM','IEUM']

    all_data = []
    for sheet in sheet_names:
        df = pd.read_excel(file_path, sheet_name=sheet)
        if y_column not in df.columns:
            print(f"Column '{y_column}' not found in sheet '{sheet}' — skipping.")
            continue
        name = label_names.get(sheet, None)
        df['Behavior'] = rf'$\texttt{{{name}}}$'  # add label to differentiate in the plot
        all_data.append(df[['maxseats', y_column, 'Behavior']])
        
    # Concatenate all data together
    if all_data:
        plot_df = pd.concat(all_data, ignore_index=True)
        
        # Define the color palette matching your dictionary
        palette = {behavior: colors[aux[behavior]] for behavior in plot_df['Behavior'].unique()}
        
        # Now plot boxplot
        ax = sns.boxplot(data=plot_df, x='maxseats', y=y_column, hue='Behavior',showfliers=False,palette=palette)
        
        num_seats = len(plot_df['maxseats'].unique())
        num_behaviors = len(plot_df['Behavior'].unique())

        for i, patch in enumerate(ax.patches):
            if i < num_behaviors*num_seats:
                behavior_index = i // num_seats
                behavior = aux_hatch[behavior_index]
                patch.set_hatch(hatch_style.get(behavior, ''))
                patch.set_edgecolor('black')
                patch.set_linewidth(2.5)
            else:
                behavior_index = (i-num_behaviors*num_seats) % num_behaviors
                behavior = aux_hatch[behavior_index]
                patch.set_hatch(hatch_style.get(behavior, ''))
                patch.set_edgecolor('black')
                patch.set_linewidth(2.5)
                
        plt.xlabel(r'Budget $B$')
        plt.ylabel(y_column_name)
        plt.xticks(ticks=range(len(all_maxseats)), labels=all_maxseats) 
        plt.ylim(bottom=0) 
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(det_stoch_path + f"boxplot_{y_column}.png", dpi=300)
        plt.close()

def plot_barplots_by_maxseats():
    """
    Plots barplots for each maxseats value, with subplots for each maxranks value.
    Each subplot compares the average of the nb of students allocated in each rank in
    the deterministic vs stochastic problem for each model (UM, CEUM, IEUM).
    """
    
    # Load data from all sheets
    data = {sheet: pd.read_excel(file_path, sheet_name=sheet) for sheet in sheet_names}

    # Get all maxseats values (assumed the same for all sheets)
    all_maxseats = sorted(data[sheet_names[0]]['maxseats'].unique())
    
    # Get all maxranks values (assumed the same for all sheets)
    all_maxranks = sorted(data[sheet_names[0]]['maxranks'].unique())

    # For each value of maxseats, create a figure with subplots for each value of maxranks
    for maxseats_val in all_maxseats:
        # Create subplots with a bit more height to accommodate the legend
        fig, axes = plt.subplots(1, len(all_maxranks), figsize=(12 * len(all_maxranks), 12), sharey=True)
        axes = np.atleast_1d(axes)  # Ensure axes is iterable 

        # Set to store unique labels for the legend
        legend_labels = set()

        for idx, maxranks_val in enumerate(all_maxranks):
            ax = axes[idx]
            x_labels = list(range(1, maxranks_val + 2))
            x = np.arange(len(x_labels))
            width = 0.25  # Bar width
            
            for i, sheet in enumerate(sheet_names):
                df = data[sheet]
                filtered = df[(df['maxseats'] == maxseats_val) & (df['maxranks'] == maxranks_val)]
                if filtered.empty:
                    continue

                # Compute average columns diff rank1 to diff rank(maxranks+1)
                avg_values = []
                for j in x_labels:
                    col_name = f'diff_rank{j}'
                    if col_name in filtered.columns:
                        avg_values.append(filtered[col_name].mean())
                    else:
                        avg_values.append(np.nan)

                # Plot: Only label bars once for each unique name
                offset = (i - 1) * width
                name = label_names.get(sheet, None)
                
                # Only add label if not already added
                if name and name not in legend_labels:
                    ax.bar(x + offset, avg_values, width=width, label=rf'$\texttt{{{name}}}$', color=colors[name], hatch=hatch_style[name], alpha=0.9, edgecolor='black',)
                    legend_labels.add(name)  # Mark label as used
                else:
                    ax.bar(x + offset, avg_values, width=width, color=colors[name], hatch=hatch_style[name], alpha=0.9, edgecolor='black',)

            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)
            ax.set_xlabel('Rank')
            ax.set_title(rf'$K = {{{maxranks_val}}}$')
            ax.grid(True)
            # ax.set_ylabel(r'$\overline{\texttt{DIFF}-\texttt{RANK}}^k$')
        axes[0].set_ylabel(r'$\overline{\texttt{DIFF-RANK}}^k$')
        fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=3)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8, bottom=0.1)

        # Save the figure
        plt.savefig(os.path.join(det_stoch_path, f'barplot_maxseats_{maxseats_val}.png'), dpi=300)
        plt.close()
        
#############################################################################################################
#                                      Comparison between behaviors                                         #
#############################################################################################################
list_comp = ['UM-CEUM','IEUM-CEUM','UM-IEUM']
# colors_comp = {'UM-CEUM':'tab:red','IEUM-CEUM':'tab:purple','UM-IEUM':'tab:brown'}
colors_comp = {'UM-CEUM': '#D55E00','IEUM-CEUM': '#CC79A7','UM-IEUM': '#56B4E9',}
style_comp = {'UM-CEUM':'-','IEUM-CEUM':'--','UM-IEUM':'-.'}
hatch_comp = {'UM-CEUM':'.','IEUM-CEUM':'-','UM-IEUM':'+'}
strat_path = save_path + 'strat/'

def plot_metric_by_budget_and_rank(type_, y_column, y_column_name):
    """
    Plots the average and quantiles of a metric (y_column) by budget, for multiple types of student behavior, and saves the figure.
    
    Parameters:
    - type_: Types: 'det' for deterministic and 'stoch' for stochastic.
    - y_column: Information to be plotted ('eval','nb_enter','nb_improv','nb_total')
    - y_column_name: Name to display on the y-axis
    """
    
    labels = {
        'UM-CEUM': ['vss_tr1', f'rp x {type_} nonstrat {y_column}(%)'],
        'IEUM-CEUM': ['vss_tr1', f'rp x {type_} simple {y_column}(%)'],
        'UM-IEUM': ['vss_tr3', f'rp x {type_} nonstrat {y_column}(%)']
    }
    
    # Load data from all sheets
    data = {label: pd.read_excel(file_path, sheet_name=info[0]) for label, info in labels.items()}
    
    # Get all maxranks values (assumed the same for all sheets)
    all_maxranks = sorted(data['UM-CEUM']['maxranks'].unique())
    
    # Get all maxseats values (assumed the same for all sheets)
    all_maxseats = sorted(data['UM-CEUM']['maxseats'].unique())

    # Initialize subplots
    fig, axes = plt.subplots(1, len(all_maxranks), figsize=(16 * len(all_maxranks), 10), sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)  # Ensure axes is iterable

    # Iterate per information
    for label, info in labels.items():
        df = data[label]
        sheet, column_name = info

        if column_name not in df.columns:
            print(f"Column '{column_name}' not found in sheet '{sheet}' — skipping.")
            continue

        agg_df = df.groupby(['maxranks', 'maxseats'], as_index=False)[column_name].mean()

        # Iterate per value of maxranks
        for rank_index, rank in enumerate(all_maxranks):
            ax = axes[rank_index]
            sub_df = agg_df[agg_df['maxranks'] == rank]

            ax.plot(sub_df['maxseats'], sub_df[column_name], marker='o',label=rf'$\texttt{{{label}}}$', color=colors_comp[label], linestyle=style_comp[label])

    for i, rank in enumerate(all_maxranks):
        ax = axes[i]
        ax.set_title(rf'$K = {{{rank}}}$')
        ax.set_xlabel(r'Budget $B$')
        ax.set_ylim(bottom=-2) 
        ax.set_ylabel(f'{y_column_name}')
        ax.set_xticks(all_maxseats)
        ax.grid(True)
        ax.legend()
    # plt.suptitle(f'Comparison of {y_column_name} across maxranks and budgets', fontsize=14)

    # Save and show
    plt.savefig(strat_path + f"comp_strat_{type_}_{y_column}.png", dpi=300)
    plt.show()
    
    #####################################################
    
    # Initialize subplots
    fig, axes = plt.subplots(1, len(all_maxranks), figsize=(16 * len(all_maxranks), 10), sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)  # Ensure axes is iterable

    # Iterate per information
    all_data = []
    for label, info in labels.items():
        df = data[label]
        sheet, column_name = info

        if column_name not in df.columns:
            print(f"Column '{column_name}' not found in sheet '{sheet}' — skipping.")
            continue

        df['Comparison'] = rf'$\texttt{{{label}}}$'  # add label to differentiate in the plot
        df = df.rename(columns={column_name: 'value'})
        all_data.append(df[['maxseats','maxranks','value','Comparison']])

    aux = {r'$\texttt{UM-CEUM}$':'UM-CEUM',r'$\texttt{UM-IEUM}$':'UM-IEUM',r'$\texttt{IEUM-CEUM}$':'IEUM-CEUM'}
    
    # Concatenate all data together
    if all_data:
        plot_df = pd.concat(all_data, ignore_index=True)
        
        # Iterate per value of maxranks
        for rank_index, rank in enumerate(all_maxranks):
            ax = axes[rank_index]
            sub_df = plot_df[plot_df['maxranks'] == rank]
            
            # Define the color palette matching your dictionary
            palette = {label: colors_comp[aux[label]] for label in sub_df['Comparison'].unique()}
            
            # Now plot boxplot
            ax_aux = sns.boxplot(data=sub_df, x='maxseats', y='value',hue='Comparison',showfliers=False,palette=palette,ax=ax)

            num_comp = len(list_comp)
            num_seats = len(all_maxseats)
            for i, patch in enumerate(ax_aux.patches):
                if i < num_comp*num_seats:
                    index = i // num_seats
                    comp = list_comp[index]
                    patch.set_hatch(hatch_comp.get(comp, ''))
                    patch.set_edgecolor('black')
                    patch.set_linewidth(2.5)
                else:
                    index = (i-num_comp*num_seats) % num_comp
                    comp = list_comp[index]
                    patch.set_hatch(hatch_comp.get(comp, ''))
                    patch.set_edgecolor('black')
                    patch.set_linewidth(2.5)
            
        for i, rank in enumerate(all_maxranks):
            ax = axes[i]
            ax.set_title(rf'$K = {{{rank}}}$')
            ax.set_xlabel(r'Budget $B$')
            ax.set_ylim(bottom=-2) 
            ax.set_ylabel(f'{y_column_name}')
            plt.xticks(ticks=range(len(all_maxseats)), labels=all_maxseats) 
            ax.grid(True)
            ax.legend()
        
        # Save and show
        plt.savefig(strat_path + f"boxplot_comp_strat_{type_}_{y_column}.png", dpi=300)
        plt.show()


#############################################################################################################
#                                      Performance different solvers                                        #
#############################################################################################################

perf_path = save_path + 'perf/'
sheet_perf = {0:'stoch_tr0',1:'stoch_tr1',3:'stoch_tr3'}

# Info. about solvers
nb_solvers = 6
solvers = [i for i in range(nb_solvers)]
optgap_solvers = [0,1,2,3]
solver_excel = ["modlconstraints","cutoffscore","asgheuristic","lagrangian _maxit300_maxitimprov15_mod0_inilp0","localsearch","simulannealing"]
solver_name  = [r'Mod. $L$-constraints',r'Cutoff score',r'$\texttt{ASG}$ heuristic',r'$\texttt{LR}$ heuristic',r'$\texttt{LS}$ heuristic',r'$\texttt{SA}$ heuristic']
# colors_perf  = ['magenta', 'dimgray', 'lightseagreen', 'midnightblue', 'darkgreen', 'darkred']
colors_perf = ['#E41A1C',  '#377EB8',  '#4DAF4A',  '#984EA3',  '#A65628',  '#999999' ]
styles_perf  = ['-', '-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))]
markers_perf = ['o', 'x', 's', '^', 'D', 'v', 'P']

def comparison_perf(small,name_excel,transform):
    """
    Plots the accumulated opt. gap, time, and ub gap by percentage of instance for the different methods 
    for the behavior {transform}, and saves the figures.
    
    Parameters:
    - name_excel: 
    """    
    
    df = pd.read_excel(name_excel, sheet_name=sheet_perf[transform])
    x = [i/df.shape[0]*100 for i in range(1,df.shape[0]+1)]
    
    ################################
    
    # Optimality gap
    plt.figure(figsize=(16, 10))
    for sol in optgap_solvers:
        y_column = 'gap1 ' + solver_excel[sol]
        if y_column not in df.columns:
            print(f"Column '{y_column}' not found in sheet '{sheet_perf[transform]}' — skipping.")
            continue
        
        y = df[y_column].sort_values()
        # plt.plot(x, y, label = solver_name[sol], color=colors_perf[sol], linestyle=styles_perf[sol], marker='o', linewidth=4)
        plt.plot(x, y, label = solver_name[sol], color=colors_perf[sol], linestyle=styles_perf[sol], linewidth=5)
    
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(perf_path + f"opt_gap_tr" + str(transform) + "_" + str(small) + ".png", dpi=300)
    plt.close() 
    
    
    ##############################
        
    # Time
    plt.figure(figsize=(16, 10))
    for sol in solvers:
        y_column = 'time_step1 ' + solver_excel[sol]
        if y_column not in df.columns:
            print(f"Column '{y_column}' not found in sheet '{sheet_perf[transform]}' — skipping.")
            continue
        
        y = df[y_column].sort_values()
        # plt.plot(x, y, label = solver_name[sol], color=colors_perf[sol], linestyle=styles_perf[sol], marker='o', linewidth=4)
        plt.plot(x, y, label = solver_name[sol], color=colors_perf[sol], linestyle=styles_perf[sol], linewidth=5)
        
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'Time (Seconds)')
    plt.legend()
    if transform == 0:
        plt.yticks([0, 1800, 3600, 5400, 7200])
    else:
        plt.yticks([0, 7200, 14400, 21600, 28800])
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(perf_path + f"time_tr" + str(transform) + "_" + str(small) + ".png", dpi=300)
    plt.close()
    
    ############################ 
    
    # Oportunity gap
    if transform != 0:
        plt.figure(figsize=(16, 10))
        for sol in solvers:
            if transform == 0 and sol == 4:
                continue
            y_column = 'gap_ub ' + solver_excel[sol]
            if y_column not in df.columns:
                print(f"Column '{y_column}' not found in sheet '{sheet_perf[transform]}' — skipping.")
                continue
            
            y = df[y_column].sort_values()
            # plt.plot(x, y, label = solver_name[sol], color=colors_perf[sol], linestyle=styles_perf[sol], marker='o', linewidth=4)
            plt.plot(x, y, label = solver_name[sol], color=colors_perf[sol], linestyle=styles_perf[sol], linewidth=5)
            
        plt.xlabel(r'Instances (\%)')
        # plt.ylabel(r'${\textnormal{gap}}_{ub^*}(\%)$')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(perf_path + f"ub_gap_tr" + str(transform) + "_" + str(small) + ".png", dpi=300)
        plt.close() 
    else:
        fig, axs = plt.subplots(2, 1, figsize=(16, 16))
        for i, exclude_bad in enumerate([False, True]):
            ax = axs[i]
            for sol in solvers:
                if exclude_bad and transform == 0 and sol == 4:
                    continue

                y_column = 'gap_ub ' + solver_excel[sol]
                if y_column not in df.columns:
                    print(f"Column '{y_column}' not found in sheet '{sheet_perf[transform]}' — skipping.")
                    continue

                y = df[y_column].sort_values()
                ax.plot(x, y,
                        label=solver_name[sol],
                        color=colors_perf[sol],
                        linestyle=styles_perf[sol],
                        linewidth=5)
            ax.grid(True)
            ax.legend()
            # ax.set_ylabel(r'${\textnormal{gap}}_{ub^*}(\%)$')
            # ax.set_xlabel(r'Instances (\%)')

        axs[1].set_xlabel(r'Instances (\%)')
        plt.tight_layout()
        plt.savefig(perf_path + f"ub_gap_tr{transform}_{small}.png", dpi=300)
        plt.close()
  
         
#############################################################################################################
#                                                Main                                                       #
#############################################################################################################

def main():
    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    os.makedirs(os.path.dirname(det_stoch_path), exist_ok=True)
    os.makedirs(os.path.dirname(strat_path), exist_ok=True)
    os.makedirs(os.path.dirname(perf_path), exist_ok=True)
    
    # Set global font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        # 'legend.title_fontsize': 44,
        'text.usetex': True,
        'axes.titlesize': 48,
        'axes.labelsize': 50,
        'xtick.labelsize': 43,
        'ytick.labelsize': 43,
        'legend.fontsize': 44,
        'font.size': 44,
        'axes.linewidth': 3,
        'lines.markersize': 10
    })
     
    # Comparison performance of the different methods
    # 0-> Nonstrategic (UM), 1-> Strategic complex (CEUM), 2-> Strategic simple (IEUM)
    for tr in [0,1,3]:  # Comparison for large instances - all behaviors
        comparison_perf('large',file_path,tr) 
    for tr in [1,3]:    # Comparison for small instances - only for strategic behaviors
        comparison_perf('small',file_path_small,tr)
    
    # Comparison stochastic vs deterministic
    plot_metric_by_budget('vss(%)', r'$\overline{\texttt{VSS}}$ (\%)')
    plot_metric_by_budget('nb_enter(%)', r'$\overline{\texttt{VSS-E}}$ (\%)')
    plot_metric_by_budget('nb_improv(%)', r'$\overline{\texttt{VSS-I}}$ (\%)')
    plot_barplots_by_maxseats()
    
    # Strategic comparison
    # Stochastic version
    plot_metric_by_budget_and_rank('stoch','eval',r"$\overline{\texttt{gap}}^{a'-a}_{N'}(\%)$")
    plot_metric_by_budget_and_rank('stoch','nb_enter',r"$\overline{\texttt{gap-e}}^{a'-a}_{N'}(\%)$")
    plot_metric_by_budget_and_rank('stoch','nb_improv',r"$\overline{\texttt{gap-i}}^{a'-a}_{N'}(\%)$")
    plot_metric_by_budget_and_rank('stoch','nb_total',r"$\overline{\texttt{gap-b}}^{a'-a}_{N'}(\%)$")
    
    # Deterministic version
    plot_metric_by_budget_and_rank('det','eval',r"$\overline{\texttt{gap}}^{a'-a}_{ev,N'}(\%)$")
    plot_metric_by_budget_and_rank('det','nb_enter',r"$\overline{\texttt{gap-e}}^{a'-a}_{ev,N'}(\%)$")
    plot_metric_by_budget_and_rank('det','nb_improv',r"$\overline{\texttt{gap-i}}^{a'-a}_{ev,N'}(\%)$")
    plot_metric_by_budget_and_rank('det','nb_total',r"$\overline{\texttt{gap-b}}^{a'-a}_{ev,N'}(\%)$")
  
        
if __name__ == "__main__":
    main()