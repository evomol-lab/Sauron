import sys

with open('static/script.js', 'r') as f:
    content = f.read()

old_rinpy_block = """
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
        }
"""

new_rinpy_block = """
        if (data.isRinpy) {
            graphViewerDiv.innerHTML = '';
            document.getElementById('sequence-container').style.display = 'none';
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
            
            // Render RinPy Tables
            if (data.rinpySummary) {
                document.getElementById('top-residues-container').style.display = 'block';
                document.getElementById('interaction-details-container').style.display = 'none'; // Keep standard network details hidden
                
                const tabsDiv = document.querySelector('.metrics-tabs');
                const contentDiv = document.querySelector('.metrics-content');
                tabsDiv.innerHTML = '';
                contentDiv.innerHTML = '';
                
                // Add Degree, Betweenness, Closeness
                const metrics = [
                    {key: 'degree', label: 'Top Degree'},
                    {key: 'betweenness', label: 'Top Betweenness'},
                    {key: 'closeness', label: 'Top Closeness'}
                ];
                
                let isFirstTab = true;
                
                metrics.forEach(m => {
                    if (data.rinpySummary[m.key] && data.rinpySummary[m.key].length > 0) {
                        const tabBtn = document.createElement('button');
                        tabBtn.className = `metric-tab ${isFirstTab ? 'active' : ''}`;
                        tabBtn.textContent = m.label;
                        tabBtn.onclick = (e) => {
                            document.querySelectorAll('.metric-tab').forEach(b => b.classList.remove('active'));
                            e.target.classList.add('active');
                            document.querySelectorAll('.metric-pane').forEach(p => p.classList.remove('active'));
                            document.getElementById(`pane-${m.key}`).classList.add('active');
                        };
                        tabsDiv.appendChild(tabBtn);
                        
                        const pane = document.createElement('div');
                        pane.className = `metric-pane ${isFirstTab ? 'active' : ''}`;
                        pane.id = `pane-${m.key}`;
                        
                        let tableHtml = '<table class="data-table"><thead><tr><th>Rank</th><th>Residue</th><th>Score</th></tr></thead><tbody>';
                        data.rinpySummary[m.key].forEach((item, index) => {
                            tableHtml += `<tr><td>${index+1}</td><td>${item.Residue}</td><td>${item.Value.toFixed(4)}</td></tr>`;
                        });
                        tableHtml += '</tbody></table>';
                        pane.innerHTML = tableHtml;
                        contentDiv.appendChild(pane);
                        
                        isFirstTab = false;
                    }
                });
                
                // Add Hinges
                if (data.rinpySummary.hinges && Object.keys(data.rinpySummary.hinges).length > 0) {
                    const tabBtn = document.createElement('button');
                    tabBtn.className = `metric-tab`;
                    tabBtn.textContent = 'Hinges (Modos Laplacianos)';
                    tabBtn.style.color = '#ff6b6b';
                    tabBtn.onclick = (e) => {
                        document.querySelectorAll('.metric-tab').forEach(b => b.classList.remove('active'));
                        e.target.classList.add('active');
                        document.querySelectorAll('.metric-pane').forEach(p => p.classList.remove('active'));
                        document.getElementById('pane-hinges').classList.add('active');
                    };
                    tabsDiv.appendChild(tabBtn);
                    
                    const pane = document.createElement('div');
                    pane.className = `metric-pane`;
                    pane.id = `pane-hinges`;
                    
                    let tableHtml = '<table class="data-table"><thead><tr><th>Mode</th><th>Hinge Residues</th></tr></thead><tbody>';
                    for (const [mode, residues] of Object.entries(data.rinpySummary.hinges)) {
                        tableHtml += `<tr><td style="font-weight:bold;color:#ff6b6b;">${mode.replace(/_/g, ' ').toUpperCase()}</td><td>${residues.join(', ')}</td></tr>`;
                    }
                    tableHtml += '</tbody></table>';
                    pane.innerHTML = tableHtml;
                    contentDiv.appendChild(pane);
                }
            }
        }
"""

content = content.replace(old_rinpy_block, new_rinpy_block)

with open('static/script.js', 'w') as f:
    f.write(content)

