import graphviz
from pathlib import Path

def create_state_diagram():
    dot = graphviz.Digraph('State_Legend', 
                          comment='Methylation and Variant States Legend',
                          graph_attr={'rankdir': 'TD'})

    # Create Legend subgraph
    with dot.subgraph(name='cluster_0') as legend:
        legend.attr(label='Legend')
        
        # Methylation states
        legend.node('A', 'Methylation States')
        legend.node('M_plus', 'M+ (Methylated)')
        legend.node('M_minus', 'M- (Unmethylated)')
        legend.edge('A', 'M_plus')
        legend.edge('A', 'M_minus')
        
        # Variant states
        legend.node('B', 'Variant States')
        legend.node('V_plus', 'V+ (Has Variant)')
        legend.node('V_minus', 'V- (No Variant)')
        legend.edge('B', 'V_plus')
        legend.edge('B', 'V_minus')
        
        # Data availability
        legend.node('C', 'Data Availability')
        legend.node('X_plus', 'X+ (Data Present)')
        legend.node('X_minus', 'X- (Data Missing)')
        legend.edge('C', 'X_plus')
        legend.edge('C', 'X_minus')

    # Create State Table subgraph
    with dot.subgraph(name='cluster_1') as states:
        states.attr(label='Variant Mapper States')
        states.node('table', '''
        State | Description
        ---|---
        M+V+ | Methylated with variant
        M-V+ | Unmethylated with variant
        M+V- | Methylated without variant
        M-V- | Unmethylated without variant
        X(M-)V+ | No methylation data, has variant
        X(V-)M+ | No variant data, methylated
        X(V-)M- | No variant data, unmethylated
        X(M-)X(V-) | No data for either
        ''', shape='plaintext')

    return dot

def save_diagram(output_dir='outputs'):
    dot = create_state_diagram()
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    dot.render(f'{output_dir}/state_legend', format='png', cleanup=True)
    dot.render(f'{output_dir}/state_legend', format='svg', cleanup=True)

if __name__ == '__main__':
    save_diagram()
