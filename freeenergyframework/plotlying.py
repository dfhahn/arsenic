import numpy as np
import seaborn as sns
import plotly.graph_objects as go

from freeenergyframework import stats


def plot_bar(df, ddg_cols, error_cols, exp_col='exp', exp_error_col='dexp', name_col='edge', title='',
             filename=None):
    '''
    Creates a plotly barplot. It takes a pandas.Dataframe df as input and plots horizontal bars grouping the values
    in the rows together. The columns which will be used are specified by ddg_cols (DDG values),
    error_cols (corresponding errors), exp_col (column with exp. values), exp_error_col (column with exp. errors)
    and name_col (column which will be used as y axis tick labels).
    '''

    # create color palette
#    colors = sns.color_palette(palette='bright')

    numEdges = df.shape[0]
    numBarsPerEdge = len(ddg_cols)
    height = 20 * (numBarsPerEdge + 0.3) * numEdges
    exp_size = height / numEdges / 2.0
    alim = np.max(np.fabs(df.loc[:, ddg_cols + [exp_col]].values) + np.fabs(
        df.loc[:, error_cols + [exp_error_col]].values)) * 1.05

    fig = go.Figure()

    # add data
    for i, (col, ecol) in enumerate(zip(ddg_cols, error_cols)):
        fig.add_trace(go.Bar(x=df.loc[:, col].values,
                             y=df[name_col].values,
                             error_x=dict(
                                 type='data',  # value of error bar given in data coordinates
                                 array=df.loc[:, ecol].values,
                                 visible=True),
                             name=col,
                             marker=dict(
#                                 color=f'rgba{colors[i]}',
                                 line=None
                             ),
                             hovertemplate = '<b>%{text}</b><extra></extra>'+
                             '<br>%{x:.2f} +- %{error_x.array:.2f} kcal mol<sup>-1</sup>',
                             text=df.loc[:,name_col].values,
                             hovertext='text',
                             orientation='h'
                             ))

    if exp_col != None:
        fig.add_trace(go.Scatter(x=df.loc[:, exp_col].values,
                                 y=df[name_col].values,
                                 name="experiment",
                                 mode='markers',
                                 marker=dict(
                                     symbol='line-ns',
                                     color='black',
                                     size=exp_size,
                                     line_width=4,
                                 ),
                                 hovertext=''
                                 ))

        fig.add_trace(go.Scatter(x=df.loc[:, exp_col].values - df.loc[:, exp_error_col].values,
                                 y=df[name_col].values,
                                 name='ExpErrors1',
                                 mode='markers',
                                 marker=dict(
                                     symbol='line-ns',
                                     color='black',
                                     size=exp_size,
                                     line_width=2,
                                 ),
                                 hovertext='',
                                 showlegend=False
                                 ))

        fig.add_trace(go.Scatter(x=df.loc[:, exp_col].values + df.loc[:, exp_error_col].values,
                                 y=df[name_col].values,
                                 name='ExpErrors2',
                                 mode='markers',
                                 hovertext='',
                                 marker=dict(
                                     symbol='line-ns',
                                     color='black',
                                     size=exp_size,
                                     line_width=2,
                                 ),
                                 showlegend=False
                                 ))

    fig.update_layout(
        title=title,
        xaxis=dict(
            title='ΔΔG [kcal mol<sup>-1</sup>]',
            titlefont_size=16,
            tickfont_size=14,
            range=(-alim, alim)
        ),
        yaxis=dict(
            title='Edge',
            titlefont_size=16,
            tickfont_size=14,
            range=(-.5, numEdges - .5)
        ),
        width=800,
        height=height,
        legend=dict(
            x=1.0,
            y=1.0,
            bgcolor='rgba(255, 255, 255, 0)',
            bordercolor='rgba(255, 255, 255, 0)',
            font_size=16
        ),
        barmode='group',
        bargap=0.3,  # gap between bars of adjacent location coordinates.
        bargroupgap=0.0  # gap between bars of the same location coordinate.
    )

    if filename is None:
        fig.show()
    elif filename.find('.html') > 0:
        fig.write_html(filename)
    else:
        fig.write_image(filename)


def _master_plot(x, y, c=None, title='',
                 xerr=None, yerr=None, names=None, smiles=None,
                 method_name_x='Experimental', method_name_y='', target_name='', plot_type='',
                 guidelines=True, origins=True,
                 statistics=['RMSE', 'MUE'], filename=None):
    nsamples = len(x)
    ax_min = min(min(x), min(y)) - 0.5
    ax_max = max(max(x), max(y)) + 0.5

    fig = go.Figure()

    # x = 0 and y = 0 axes through origin
    if origins:
        # x=0
        fig.add_trace(go.Scatter(x=[0, 0],
                                 y=[ax_min, ax_max],
                                 line_color='grey',
                                 mode='lines',
                                 showlegend=False
                                 ))
        # y =0
        fig.add_trace(go.Scatter(x=[ax_min, ax_max],
                                 y=[0, 0],
                                 line_color='grey',
                                 mode='lines',
                                 showlegend=False
                                 ))
    if guidelines:
        small_dist = 0.5
        fig.add_trace(go.Scatter(x=[ax_min, ax_max, ax_max, ax_min],
                                 y=[ax_min + 2. * small_dist, ax_max + 2. * small_dist, ax_max - 2. * small_dist,
                                    ax_min - 2. * small_dist],
                                 name='1 kcal/mol margin',
                                 hoveron='points+fills',
                                 hoverinfo='name',
                                 fill='toself',
                                 mode='lines', line_width=0, fillcolor='rgba(0, 0, 0, 0.2)',
                                 showlegend=False))

        fig.add_trace(go.Scatter(x=[ax_min, ax_max, ax_max, ax_min],
                                 y=[ax_min + small_dist, ax_max + small_dist, ax_max - small_dist, ax_min - small_dist],
                                 name='.5 kcal/mol margin',
                                 hoveron='points+fills',
                                 hoverinfo='name',
                                 fill='toself',
                                 mode='lines', line_width=0, fillcolor='rgba(0, 0, 0, 0.2)',
                                 showlegend=False))

    # diagonal
    fig.add_trace(go.Scatter(x=[ax_min, ax_max],
                             y=[ax_min, ax_max],
                             line_color='black',
                             line_dash='dash',
                             mode='lines',
                             showlegend=False
                             ))

    # 2.372 kcal / mol = 4 RT
    if c is not None:
        clr = c
        cbar =  dict(
                                     title="Analytical Error<br>[kcal mol<sup>-1</sup>]",
                                     # lenmode="pixels", len=200,
                                     thicknessmode='pixels', thickness=10,
                                     yanchor="top", y=1,
                                     tickvals=[0, 5, 7.5],
                                     ticktext=['0.0', '5.0', '>5.0']
                                 )
    else:
        clr = np.zeros_like(x)
        cbar = None

    if names is not None:
        if smiles is not None:
            smile_strings = ['<br>'.join(smile) for smile in smiles]
            n = [f'{name}<br>{smile}' for name, smile in zip(names, smile_strings)]
        else:
            n = names
    else:
        n = ['' for x in clr]

    fig.add_trace(go.Scatter(x=x, y=y,
                             mode='markers',
                             name=f'{target_name},{method_name_y}',
                             marker=dict(
                                 symbol='circle',
                                 color=clr,
                                 colorscale=[(0, 'rgb(0,0,255)'), (0.5, 'rgb(127,0, 127)'), (0.5, 'rgb(255, 0, 0)'),
                                             (1.0, 'rgb(255,0,0)')],
                                 cmax=10.0,
                                 cmin=0.0,
                                 colorbar=cbar
                             ),
                             error_x=dict(
                                 type='data',  # value of error bar given in data coordinates
                                 array=xerr,
                                 visible=True),
                             error_y=dict(
                                 type='data',  # value of error bar given in data coordinates
                                 array=yerr,
                                 visible=True),
                             hovertemplate='<b>%{text}</b><extra></extra>' +
                                           '<br>Exp: %{x:.2f} +- %{error_x.array:.2f} kcal mol<sup>-1</sup>' +
                                           '<br>Calc: %{y:.2f} +- %{error_y.array:.2f} kcal mol<sup>-1</sup>',
                             text=n,
                             hovertext='text',
                             showlegend=False
                             ))

    # stats and title
    string = []
    for statistic in statistics:
        s = stats.bootstrap_statistic(x, y, statistic=statistic)
        string.append(
            f"{statistic + ':':5s}{s['mle']:5.2f} [95%: {s['low']:5.2f}, {s['high']:5.2f}] kcal mol<sup>-1</sup>")
    statistics_string = '<br>'.join(string)

    long_title = f'{title}<br>{target_name} (N = {nsamples})<br>{statistics_string}'

    # figure layout
    fig.update_layout(
        title=dict(
            text=long_title,
            font_family='monospace',
            x=0.0,
            y=0.99,
            font_size=14,
        ),
        xaxis=dict(
            title=f'{method_name_x} {plot_type} [kcal mol<sup>-1</sup>]',
            titlefont_size=14,
            tickfont_size=12,
            range=(ax_min, ax_max)
        ),
        yaxis=dict(
            title=f'{method_name_y} {plot_type}  [kcal mol<sup>-1</sup>]',
            titlefont_size=14,
            tickfont_size=12,
            range=(ax_min, ax_max)
        ),
        width=400,
        height=400
        #         legend=dict(
        #             x=1.0,
        #             y=1.0,
        #             bgcolor='rgba(255, 255, 255, 0)',
        #             bordercolor='rgba(255, 255, 255, 0)',
        #             font_size=12
        #         )
    )

    if not np.any(c):
        fig.update_layout()

    if filename is None:
        fig.show()
    elif filename.find('.html') > 0:
        fig.write_html(filename)
    else:
        fig.write_image(filename)
