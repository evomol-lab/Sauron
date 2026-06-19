import sys

with open('static/script.js', 'r') as f:
    content = f.read()

old_graph_load = """
        // Load Graph
        if (data.edgesFile && data.nodesFile) {
"""

new_graph_load = """
        // Load Graph
        if (data.isRinpy) {
            graphViewerDiv.innerHTML = '';
            document.getElementById('interaction-details-container').style.display = 'none';
            document.getElementById('sequence-container').style.display = 'none';
            document.getElementById('top-residues-container').style.display = 'none';
            document.getElementById('interchain-container').style.display = 'none';
            document.getElementById('ligand-container').style.display = 'none';
            
            if (data.rinpyHtml) {
                const iframe = document.createElement('iframe');
                iframe.src = `/rinpy_view/${data.rinpyPrefix}/${data.rinpyHtml}`;
                iframe.style.width = '100%';
                iframe.style.height = '100%';
                iframe.style.border = 'none';
                iframe.style.borderRadius = '8px';
                graphViewerDiv.appendChild(iframe);
            } else {
                graphViewerDiv.innerHTML = '<div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--text-muted);text-align:center;">RinPy calculation completed.<br><br>Download the ZIP file to view the full CSV reports and interactive plots.</div>';
            }
        } else if (data.edgesFile && data.nodesFile) {
"""

content = content.replace(old_graph_load, new_graph_load)

with open('static/script.js', 'w') as f:
    f.write(content)

