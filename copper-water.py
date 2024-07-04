from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
import dash

app = dash.Dash(__name__)

app.layout = html.Div([
    html.Div([
        html.H1(u'Copper-H\u2082O system')],
        style={
            'text-align':'center',
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '-4px 0 -4px 0',
            'margin-bottom': '2px',
            'color': '#10328f',
        }
    ),

    html.Div([
        html.H6(u"Total copper(I) concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='copper_slider1',
            min=0.1,
            max=2.0,
            value=1.0,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]},
        ),
    ],
        style={
            'padding': '10px 15px 10px',
            'border': 'thin lightgrey solid',
            'margin-bottom': '3px',
            'backgroundColor': 'rgb(250, 250, 250)',
        }
    ),

    html.Div([
        html.H6(u"Total copper(II) concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='copper_slider2',
            min=0.1,
            max=2.0,
            value=1.0,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]},
        ),
    ],
        style={
            'padding': '10px 15px 10px',
            'border': 'thin lightgrey solid',
            'margin-bottom': '3px',
            'backgroundColor': 'rgb(250, 250, 250)',
        }
    ),

    html.Div([
        html.Div([
            dcc.Graph(id='copperpure1_speciation')
        ],
        style={
            'width': '48%',
            'margin-left': '1%',
        }
        ),
        html.Div([
            dcc.Graph(id='copperpure_potential')
        ],
        style={
            'width': '48%',
            'margin-right': '1%',
        }
        )
    ],
    className='row',
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '0 0px 0 15px',
            "margin-left": "0px",
            'margin-right': '0px',
            'margin-bottom': '3px'
        }),
    html.Div([
        html.Div([
            dcc.Graph(id='copperpure2_speciation')
        ],
            style={
                'width': '48%',
                'margin-left': '1%',
            }
        ),
        html.Div([
            dcc.Graph(id='copperpure3_speciation')
        ],
            style={
                'width': '48%',
                'margin-right': '1%',
            }
        )
    ],
        className='row',
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '0 0px 0 15px',
            "margin-left": "0px",
            'margin-right': '0px',
            'margin-bottom': '3px'
        }),
    ],

    style = {
        'margin-top': '40px',
        'margin-bottom': '5px'
    }
)

@app.callback(
    Output('copperpure_potential', 'figure'),
    [Input('copper_slider2', 'value')])
def speciation_graph(Cu_total):
    cu2p = cuo2 = Cu_total
    T_ = 298
    pH_x = np.linspace(-2, 16, 71)

    def trace_generator(pH_x, cu2p, cuo2, T_):

        def interceptgenerator(pH_x, cu2p, cuo2, T_):

            bottom_eqs = [[R1(cu2p, T_) for i in range(len(pH_x))], R2(pH_x, T_, cu2p), R3(pH_x, T_), R4(pH_x, T_, cuo2)]
            y_vars = []
            for i, eq in enumerate(bottom_eqs):
                y_vars.append(np.polyfit(pH_x, eq, 1))

            As = []
            Bs = []
            for i, var in enumerate(y_vars):
                if i >= 1:
                    As.append(np.array([[-y_vars[i - 1][0], 1], [-y_vars[i][0], 1]]))
                    Bs.append(np.array([y_vars[i - 1][1], y_vars[i][1]]))
            inters = []
            for i, ms in enumerate(As):
                for j, cs in enumerate(Bs):
                    if i == j:
                        inters.append(np.linalg.inv(As[i]).dot(Bs[j]))
            return inters

        inters = interceptgenerator(pH_x, cu2p, cuo2, T_)
        xinters = []
        for item in inters:
            xinters.append(item[0])

        x_data = []
        for i, item in enumerate(xinters):
            if i == 0:
                x_data.append(np.linspace(-2, item, 5))
            elif i >= 1:
                x_data.append(np.linspace(xinters[i - 1], item, 5))
        finalindex = len(xinters) - 1
        x_data.append(np.linspace(xinters[finalindex], 16, 5))

        y_data_bottom = [[R1(cu2p, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_, cu2p), R3(x_data[2], T_), R4(x_data[3], T_, cuo2)]
        new_x_bottom = []
        new_y_bottom = []

        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_bottom.append(xvalue)
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        x_data_verticals_1 = np.linspace(inters[1][0], inters[1][0], 5)
        y_data_verticals_1 = np.linspace(inters[1][1], 2.4, 5)
        x_data_verticals_2 = np.linspace(inters[2][0], inters[2][0], 5)
        y_data_verticals_2 = np.linspace(inters[2][1], 2.4, 5)
        T1_xdata = np.linspace(inters[0][0], 16, 5)
        T1_ydata = np.linspace(inters[0][1], B1(16, 298), 5)
        abx = np.linspace(0, x_data[0], 5)

        cu2pregionx = list(x_data[0]) + list(x_data[1]) + list(x_data_verticals_1) + list(np.linspace(-2, inters[1][0], 5)) + list([-2 for i in range(0, 5)])
        cu2pregiony = list(y_data_bottom[0]) + list(y_data_bottom[1]) + list(y_data_verticals_1) + list([2.4 for i in range(0, 5)]) + list(reversed(np.linspace(inters[0][1], 2.4, 5)))

        cuoregionx = list(x_data_verticals_1) + list(x_data[2]) + list(x_data_verticals_2) + list(reversed(x_data[2]))
        cuoregiony = list(y_data_verticals_1) + list(y_data_bottom[2]) + list(y_data_verticals_2) + list([2.4 for i in range(0, 5)])

        cuo2regionx = list(x_data_verticals_2) + list(x_data[3]) + list([16 for i in range(0, 5)]) + list(reversed(x_data[3]))
        cuo2regiony = list(y_data_verticals_2) + list(y_data_bottom[3]) + list(np.linspace(inters[2][1], 2.4, 5)) + list([2.4 for i in range(0, 5)])

        trace0 = go.Scatter(
            x=new_x_bottom,
            y=new_y_bottom,
            mode='lines',
            name="x",
            line_color='#2E86AB'
        )

        trace1 = go.Scatter(
            x=x_data_verticals_1,
            y=y_data_verticals_1,
            mode='lines',
            name="y",
            line_color='#2E86AB'
        )

        trace2 = go.Scatter(
            x=x_data_verticals_2,
            y=y_data_verticals_2,
            mode='lines',
            name="y",
            line_color='#2E86AB'
        )

        trace3 = go.Scatter(
            x=T1_xdata,
            y=T1_ydata,
            mode='lines',
            name="y",
            line_color='#2E86AB'
        )

        cu2pregions = go.Scatter(
            x=cu2pregionx,
            y=cu2pregiony,
            fill='toself',
            mode='none',
            fillcolor='#f6c3be',
            showlegend=False,
            hoveron='fills',
            text='Cu<sub>2</sub>O',
            hoverinfo='text'
        )

        cuoregions = go.Scatter(
            x=cuoregionx,
            y=cuoregiony,
            fill='toself',
            mode='none',
            fillcolor='#ffebbb',
            showlegend=False,
            hoveron='fills',
            text='Cu<sup>o</sup>',
            hoverinfo='text'
        )

        cuo2regions = go.Scatter(
            x=cuo2regionx,
            y=cuo2regiony,
            fill='toself',
            mode='none',
            fillcolor='#c8d4e3',
            showlegend=False,
            hoveron='fills',
            text='Cu(OH)<sub>2</sub>',
            hoverinfo='text'
        )

        data = [cu2pregions, cuoregions, cuo2regions, trace0, trace1, trace2, trace3]
        return data

    data = trace_generator(pH_x, cu2p, cuo2, T_)
    layout = go.Layout(
        showlegend=False,
        xaxis=dict(
            showline=True,
            linecolor='black',
            ticks='inside',
            showgrid=True,
            mirror=True,
            gridcolor='#E6E6E6',
            title=r'$ \text{pH}$',
            range=[-2, 16],
            zeroline=False
        ),
        yaxis=dict(
            showline=True,
            linecolor='black',
            ticks='inside',
            showgrid=True,
            mirror=True,
            gridcolor='#E6E6E6',
            title=r'$E  \text{ (V) vs. SCE}$',
            range=[-1, 2.4],
            zeroline=False
        ),
        margin=dict(t=20, b=10, r=10),
        width=450,
        height=450
    )

    fig = go.Figure(data=data, layout=layout)
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
