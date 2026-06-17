document.addEventListener('DOMContentLoaded', () => {
    const themeToggle = document.getElementById('theme-toggle');
    const logoImg = document.getElementById('logo-img');
    const themeIconLight = document.getElementById('theme-icon-light');
    const themeIconDark = document.getElementById('theme-icon-dark');
    
    // Check saved theme
    if (localStorage.getItem('theme') === 'light') {
        document.body.classList.add('light-mode');
        logoImg.src = '/static/sauron-rins-logo.png';
        themeIconDark.style.display = 'none';
        themeIconLight.style.display = 'block';
    }

    themeToggle.addEventListener('click', () => {
        document.body.classList.toggle('light-mode');
        if (document.body.classList.contains('light-mode')) {
            localStorage.setItem('theme', 'light');
            logoImg.src = '/static/sauron-rins-logo.png';
            themeIconDark.style.display = 'none';
            themeIconLight.style.display = 'block';
        } else {
            localStorage.setItem('theme', 'dark');
            logoImg.src = '/static/sauron-rins-logo-inv.png';
            themeIconLight.style.display = 'none';
            themeIconDark.style.display = 'block';
        }
    });

    const dropZone = document.getElementById('drop-zone');
    const fileInput = document.getElementById('file-input');
    const fileNameDisplay = document.getElementById('file-name');
    const submitBtn = document.getElementById('submit-btn');
    const uploadForm = document.getElementById('upload-form');
    const errorMessage = document.getElementById('error-message');
    const successState = document.getElementById('success-state');
    const downloadBtn = document.getElementById('download-btn');
    const resetBtn = document.getElementById('reset-btn');

    let selectedFile = null;
    
    window.toggleInputMode = () => {
        const mode = document.querySelector('input[name="input_mode"]:checked').value;
        const dropZone = document.getElementById('drop-zone');
        const fetchZone = document.getElementById('fetch-zone');
        const submitBtn = document.getElementById('submit-btn');
        
        if (mode === 'upload') {
            dropZone.style.display = 'flex';
            fetchZone.style.display = 'none';
            submitBtn.disabled = !selectedFile;
        } else {
            dropZone.style.display = 'none';
            fetchZone.style.display = 'flex';
            window.validateFetchInput();
        }
    };
    
    window.validateFetchInput = () => {
        const mode = document.querySelector('input[name="input_mode"]:checked').value;
        if (mode === 'fetch') {
            const fetchId = document.getElementById('fetch_id').value.trim();
            document.getElementById('submit-btn').disabled = fetchId === '';
        }
    };

    // Handle Drag & Drop
    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
        dropZone.addEventListener(eventName, preventDefaults, false);
    });

    function preventDefaults(e) {
        e.preventDefault();
        e.stopPropagation();
    }

    ['dragenter', 'dragover'].forEach(eventName => {
        dropZone.addEventListener(eventName, () => dropZone.classList.add('dragover'), false);
    });

    ['dragleave', 'drop'].forEach(eventName => {
        dropZone.addEventListener(eventName, () => dropZone.classList.remove('dragover'), false);
    });

    dropZone.addEventListener('drop', (e) => {
        const dt = e.dataTransfer;
        const files = dt.files;
        handleFiles(files);
    });

    // Handle Click Upload
    dropZone.addEventListener('click', () => fileInput.click());
    
    fileInput.addEventListener('change', function() {
        handleFiles(this.files);
    });

    function handleFiles(files) {
        if (files.length > 0) {
            const file = files[0];
            const validExtensions = ['.pdb', '.cif'];
            const fileExt = file.name.substring(file.name.lastIndexOf('.')).toLowerCase();
            
            if (validExtensions.includes(fileExt)) {
                selectedFile = file;
                fileNameDisplay.textContent = file.name;
                submitBtn.disabled = false;
                errorMessage.style.display = 'none';
            } else {
                fileNameDisplay.textContent = 'Invalid file format. Please upload .pdb or .cif';
                selectedFile = null;
                submitBtn.disabled = true;
            }
        }
    }

    // Form Submission via Fetch API
    uploadForm.addEventListener('submit', async (e) => {
        e.preventDefault();
        
        const mode = document.querySelector('input[name="input_mode"]:checked').value;
        if (mode === 'upload' && !selectedFile) return;
        if (mode === 'fetch' && document.getElementById('fetch_id').value.trim() === '') return;

        // UI State: Loading
        submitBtn.classList.add('loading');
        submitBtn.disabled = true;
        errorMessage.style.display = 'none';

        const formData = new FormData();
        formData.append('input_mode', mode);
        
        if (mode === 'upload') {
            formData.append('file', selectedFile);
        } else {
            formData.append('fetch_db', document.getElementById('fetch_db').value);
            formData.append('fetch_id', document.getElementById('fetch_id').value.trim());
        }
        formData.append('calc_method', document.getElementById('calc_method').value);
        formData.append('strict_angle', document.getElementById('strict_angle').checked);
        formData.append('remove_multiples', document.getElementById('remove_multiples').checked);
        
        const chainsOpt = document.querySelector('input[name="chains_opt"]:checked').value;
        formData.append('chains_opt', chainsOpt);
        formData.append('chains_input', document.getElementById('chains_input').value);

        try {
            const response = await fetch('/upload', {
                method: 'POST',
                body: formData
            });

            const data = await response.json();

            if (!response.ok) {
                throw new Error(data.error || 'Computation failed.');
            }

            // UI State: Success
            uploadForm.style.display = 'none';
            successState.style.display = 'block';
            document.querySelector('.container').classList.add('expanded');
            downloadBtn.href = data.downloadUrl;

            // Load Visualizations
            loadVisualizations(data);

        } catch (error) {
            // UI State: Error
            errorMessage.textContent = error.message;
            errorMessage.style.display = 'block';
        } finally {
            submitBtn.classList.remove('loading');
            submitBtn.disabled = false;
        }
    });

    // Reset Form
    resetBtn.addEventListener('click', () => {
        successState.style.display = 'none';
        uploadForm.style.display = 'block';
        document.querySelector('.container').classList.remove('expanded');
        document.getElementById('structure-viewer').innerHTML = '';
        document.getElementById('graph-viewer').innerHTML = '';
        document.getElementById('interaction-details-container').style.display = 'none';
        document.getElementById('sequence-container').style.display = 'none';
        document.getElementById('sequence-content').innerHTML = '';
        document.getElementById('top-residues-container').style.display = 'none';
        document.getElementById('interchain-container').style.display = 'none';
        document.getElementById('ligand-container').style.display = 'none';
        document.getElementById('complexity-analysis-container').style.display = 'none';
        document.querySelector('#interaction-table tbody').innerHTML = '';
        document.querySelector('#top-residues-table tbody').innerHTML = '';
        document.querySelector('#interchain-table tbody').innerHTML = '';
        document.querySelector('#ligand-table tbody').innerHTML = '';
        
        const existingWarning = document.getElementById('multimodel-warning');
        if (existingWarning) existingWarning.remove();
        
        uploadForm.reset();
        selectedFile = null;
        fileNameDisplay.textContent = 'No file selected';
        submitBtn.disabled = true;
        errorMessage.style.display = 'none';
    });

    async function loadVisualizations(data) {
        const structViewerDiv = document.getElementById('structure-viewer');
        const graphViewerDiv = document.getElementById('graph-viewer');
        
        const existingWarning = document.getElementById('multimodel-warning');
        if (existingWarning) existingWarning.remove();
        
        if (data.isMultimodel) {
            const warningDiv = document.createElement('div');
            warningDiv.id = 'multimodel-warning';
            warningDiv.style.padding = '12px 16px';
            warningDiv.style.backgroundColor = 'rgba(239, 68, 68, 0.1)';
            warningDiv.style.border = '1px solid var(--error-color)';
            warningDiv.style.borderRadius = '8px';
            warningDiv.style.marginBottom = '20px';
            warningDiv.style.color = 'var(--text-main)';
            warningDiv.style.lineHeight = '1.5';
            warningDiv.innerHTML = `<strong>Multimodel structure file detected.</strong> Only the visualization features for the first model will be shown. For comparing RINs across conformations, use our other tool <a href="https://github.com/evomol-lab/slytherins" target="_blank" style="color: var(--primary-color); font-weight: 600;">SlyTheRINs</a>.`;
            document.querySelector('.success-state').insertBefore(warningDiv, document.querySelector('.visualization-panel'));
        }
        
        structViewerDiv.innerHTML = '<div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--text-muted)">Loading structure...</div>';
        graphViewerDiv.innerHTML = '<div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--text-muted)">Loading graph...</div>';

        let molstarViewer = null;
        let graph = null;

        // Load 3D Structure (PDBe Molstar)
        try {
            structViewerDiv.innerHTML = '';
            
            // Options for PDBe Molstar
            const options = {
                customData: {
                    url: `/view/${data.pdbFile}`,
                    format: data.pdbFile.endsWith('.cif') ? 'cif' : 'pdb'
                },
                hideControls: true,
                bgColor: {r: 0, g: 0, b: 0},
                lighting: 'plastic',
                visualStyle: 'cartoon'
            };
            
            molstarViewer = new PDBeMolstarPlugin();
            molstarViewer.render(structViewerDiv, options);
            
        } catch (e) {
            console.error("Failed to load structure", e);
            structViewerDiv.innerHTML = '<div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--error-color)">Failed to load structure</div>';
        }

        // Load Graph
        if (data.edgesFile && data.nodesFile) {
            try {
                const requests = [
                    fetch(`/view/${data.edgesFile}`),
                    fetch(`/view/${data.nodesFile}`)
                ];
                if (data.metricsFile) {
                    requests.push(fetch(`/view/${data.metricsFile}`));
                }
                
                const responses = await Promise.all(requests);
                const edgesText = await responses[0].text();
                const nodesText = await responses[1].text();
                const metricsText = data.metricsFile ? await responses[2].text() : '';

                const graphData = parseNetworkData(nodesText, edgesText);
                const metricsData = parseMetricsData(metricsText);
                
                // Inject degree into nodes for sizing
                graphData.nodes.forEach(node => {
                    const m = metricsData[node.id];
                    node.degree = m && m.degree !== 'N/A' ? parseInt(m.degree) : 1;
                });
                
                const interactionColors = {
                    'HBOND': '#3b82f6', // blue
                    'VDW': '#9ca3af',   // gray
                    'IONIC': '#ef4444', // red
                    'PIPISTACK': '#10b981', // green
                    'PICATION': '#f59e0b', // orange
                    'SSBOND': '#eab308' // yellow
                };
                
                const dsspColors = {
                    'H': '#e85d75', // Alpha helix
                    'B': '#ffb347', // Beta bridge
                    'E': '#f4d35e', // Strand
                    'G': '#f08a5d', // 3-10 helix
                    'I': '#db5461', // Pi helix
                    'T': '#4db6ac', // Turn
                    'S': '#81c784', // Bend
                    'C': '#e0e0e0', // Coil
                    ' ': '#e0e0e0'
                };

                const chemicalColors = {
                    'ALA': '#9ca3af', 'VAL': '#9ca3af', 'LEU': '#9ca3af', 'ILE': '#9ca3af', 'MET': '#9ca3af', 'TRP': '#9ca3af', 'PHE': '#9ca3af', 'PRO': '#9ca3af', // Non-polar
                    'SER': '#10b981', 'THR': '#10b981', 'CYS': '#10b981', 'TYR': '#10b981', 'ASN': '#10b981', 'GLN': '#10b981', // Polar
                    'LYS': '#3b82f6', 'ARG': '#3b82f6', 'HIS': '#3b82f6', // Basic
                    'ASP': '#ef4444', 'GLU': '#ef4444' // Acidic
                };
                
                const chainColors = {};
                const palette = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'];
                graphData.nodes.forEach(node => {
                    if (node.chain && !chainColors[node.chain]) {
                        chainColors[node.chain] = palette[Object.keys(chainColors).length % palette.length];
                    }
                });

                const getNodeColor = (node, type) => {
                    if (type === 'dssp') return dsspColors[node.dssp] || '#e0e0e0';
                    if (type === 'chemical') return chemicalColors[node.residue] || '#e0e0e0';
                    if (type === 'chain') return chainColors[node.chain] || '#e0e0e0';
                    return '#6366f1';
                };
                
                // Edge Legend
                const edgeLegendContainer = document.getElementById('edge-legend');
                edgeLegendContainer.innerHTML = '';
                for (const key in interactionColors) {
                    const item = document.createElement('div');
                    item.style.display = 'flex';
                    item.style.alignItems = 'center';
                    item.style.gap = '4px';
                    item.innerHTML = `<span style="display:inline-block; width:12px; height:3px; background-color:${interactionColors[key]};"></span><span>${key}</span>`;
                    edgeLegendContainer.appendChild(item);
                }
                
                const updateLegend = (type) => {
                    const legendContainer = document.getElementById('color-legend');
                    legendContainer.innerHTML = '';
                    
                    let colors = {};
                    let labels = {};
                    
                    if (type === 'dssp') {
                        colors = dsspColors;
                        labels = {
                            'H': 'Alpha helix', 'B': 'Beta bridge', 'E': 'Strand', 
                            'G': '3-10 helix', 'I': 'Pi helix', 'T': 'Turn', 
                            'S': 'Bend', 'C': 'Coil', ' ': 'Unknown'
                        };
                    } else if (type === 'chemical') {
                        colors = {
                            'Non-polar': '#9ca3af', 'Polar': '#10b981', 
                            'Basic': '#3b82f6', 'Acidic': '#ef4444'
                        };
                        labels = {
                            'Non-polar': 'Non-polar', 'Polar': 'Polar', 
                            'Basic': 'Basic', 'Acidic': 'Acidic'
                        };
                    } else if (type === 'chain') {
                        colors = chainColors;
                        for (let c in chainColors) {
                            labels[c] = `Chain ${c}`;
                        }
                    } else {
                        return;
                    }
                    
                    for (const key in colors) {
                        const color = colors[key];
                        const label = labels[key];
                        if (label && (color !== '#e0e0e0' || key === 'C')) {
                            const item = document.createElement('div');
                            item.style.display = 'flex';
                            item.style.alignItems = 'center';
                            item.style.gap = '4px';
                            item.innerHTML = `<span style="display:inline-block; width:10px; height:10px; border-radius:50%; background-color:${color};"></span><span>${label}</span>`;
                            legendContainer.appendChild(item);
                        }
                    }
                };
                
                graphViewerDiv.innerHTML = '';
                const initColorType = document.getElementById('node-color-select').value;
                updateLegend(initColorType);
                
                graph = ForceGraph3D()(graphViewerDiv)
                    .graphData(graphData)
                    .nodeLabel('id')
                    .nodeVal(node => node.degree || 1)
                    .nodeColor(node => getNodeColor(node, initColorType))
                    .linkColor(link => {
                        if (link.interaction) {
                            const type = link.interaction.split(':')[0];
                            return interactionColors[type] || '#9ca3af';
                        }
                        return '#9ca3af';
                    })
                    .linkWidth(2)
                    .backgroundColor('rgba(0,0,0,0)')
                    .onNodeClick(node => {
                        // Focus camera on node
                        const distance = 40;
                        const distRatio = 1 + distance/Math.hypot(node.x, node.y, node.z);
                        graph.cameraPosition(
                            { x: node.x * distRatio, y: node.y * distRatio, z: node.z * distRatio },
                            node,
                            1000 // ms transition duration
                        );
                        
                        // Highlight in Molstar
                        if (molstarViewer) {
                            molstarViewer.visual.clearSelection();
                            molstarViewer.visual.select({
                                data: [{
                                    struct_asym_id: node.chain,
                                    start_residue_number: parseInt(node.position),
                                    end_residue_number: parseInt(node.position)
                                }]
                            });
                            
                            molstarViewer.visual.focus([{
                                struct_asym_id: node.chain,
                                start_residue_number: parseInt(node.position),
                                end_residue_number: parseInt(node.position)
                            }]);
                        }
                        
                        // Update Interaction Table
                        updateInteractionTable(node, graphData.links, metricsData);
                    });
                    
                document.getElementById('node-color-select').addEventListener('change', (e) => {
                    const sel = e.target.value;
                    graph.nodeColor(node => getNodeColor(node, sel));
                    updateLegend(sel);
                });
                
                // Top 10 Residues Logic
                const renderTopResidues = () => {
                    const container = document.getElementById('top-residues-container');
                    if (Object.keys(metricsData).length === 0) {
                        container.style.display = 'none';
                        return;
                    }
                    
                    container.style.display = 'block';
                    const select = document.getElementById('top-feature-select');
                    const feature = select.value;
                    const tbody = document.querySelector('#top-residues-table tbody');
                    const featureHeader = document.getElementById('top-feature-header');
                    
                    const featureLabels = {
                        'degree': 'Degree',
                        'clustering': 'Clustering Coefficient',
                        'betweenness': 'Betweenness Centrality',
                        'eigenvector': 'Eigenvector Centrality'
                    };
                    featureHeader.textContent = featureLabels[feature];
                    
                    const nodesArray = Object.keys(metricsData).map(id => {
                        return {
                            id: id,
                            value: metricsData[id][feature] !== 'N/A' ? parseFloat(metricsData[id][feature]) : NaN
                        };
                    }).filter(n => !isNaN(n.value));
                    
                    nodesArray.sort((a, b) => b.value - a.value);
                    
                    const top10 = nodesArray.slice(0, 10);
                    
                    tbody.innerHTML = '';
                    top10.forEach(n => {
                        const tr = document.createElement('tr');
                        tr.style.borderBottom = '1px solid var(--panel-border)';
                        tr.style.cursor = 'pointer';
                        
                        tr.addEventListener('mouseenter', () => {
                            tr.style.backgroundColor = 'rgba(255, 255, 255, 0.05)';
                        });
                        tr.addEventListener('mouseleave', () => {
                            tr.style.backgroundColor = 'transparent';
                        });
                        
                        tr.innerHTML = `
                            <td style="padding: 8px; color: var(--primary-color);">${n.id}</td>
                            <td style="padding: 8px;">${feature === 'degree' ? n.value : n.value.toFixed(4)}</td>
                        `;
                        
                        tr.addEventListener('click', () => {
                            const graphNode = graphData.nodes.find(gn => gn.id === n.id);
                            if (graphNode) {
                                // Simulate focus camera
                                const distance = 40;
                                const distRatio = 1 + distance/Math.hypot(graphNode.x, graphNode.y, graphNode.z);
                                graph.cameraPosition(
                                    { x: graphNode.x * distRatio, y: graphNode.y * distRatio, z: graphNode.z * distRatio },
                                    graphNode,
                                    1000
                                );
                                
                                // Highlight in Molstar
                                if (molstarViewer) {
                                    molstarViewer.visual.clearSelection();
                                    molstarViewer.visual.select({
                                        data: [{
                                            struct_asym_id: graphNode.chain,
                                            start_residue_number: parseInt(graphNode.position),
                                            end_residue_number: parseInt(graphNode.position)
                                        }]
                                    });
                                    
                                    molstarViewer.visual.focus([{
                                        struct_asym_id: graphNode.chain,
                                        start_residue_number: parseInt(graphNode.position),
                                        end_residue_number: parseInt(graphNode.position)
                                    }]);
                                }
                                
                                updateInteractionTable(graphNode, graphData.links, metricsData);
                            }
                        });
                        
                        tbody.appendChild(tr);
                    });
                };
                
                renderTopResidues();
                document.getElementById('top-feature-select').addEventListener('change', renderTopResidues);
                
                // Primary Sequence Representation
                const renderSequence = () => {
                    const seqContainer = document.getElementById('sequence-container');
                    const seqContent = document.getElementById('sequence-content');
                    
                    if (!graphData.nodes || graphData.nodes.length === 0) {
                        return;
                    }
                    
                    const aminoAcidMap = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                    };
                    
                    const dsspLabels = {
                        'H': 'Alpha helix', 'B': 'Beta bridge', 'E': 'Strand', 
                        'G': '3-10 helix', 'I': 'Pi helix', 'T': 'Turn', 
                        'S': 'Bend', 'C': 'Coil', ' ': 'Unknown'
                    };
                    
                    const seqLegend = document.getElementById('sequence-legend');
                    seqLegend.innerHTML = '<span style="font-weight: 600; margin-right: 8px;">DSSP Structure:</span>';
                    for (const k in dsspLabels) {
                        const color = dsspColors[k] || '#e0e0e0';
                        const item = document.createElement('div');
                        item.style.display = 'flex';
                        item.style.alignItems = 'center';
                        item.style.gap = '4px';
                        item.innerHTML = `<span style="display:inline-block; width:12px; height:12px; border-radius:50%; background-color:${color};"></span><span>${dsspLabels[k]}</span>`;
                        seqLegend.appendChild(item);
                    }
                    
                    const chains = {};
                    graphData.nodes.forEach(node => {
                        if (!chains[node.chain]) {
                            chains[node.chain] = [];
                        }
                        chains[node.chain].push(node);
                    });
                    
                    seqContent.innerHTML = '';
                    
                    Object.keys(chains).sort().forEach(chain => {
                        const nodesInChain = chains[chain].sort((a, b) => parseInt(a.position) - parseInt(b.position));
                        
                        const chainDiv = document.createElement('div');
                        chainDiv.style.display = 'flex';
                        chainDiv.style.flexDirection = 'column';
                        chainDiv.style.gap = '8px';
                        
                        const chainLabel = document.createElement('div');
                        chainLabel.textContent = `Chain ${chain}`;
                        chainLabel.style.fontWeight = 'bold';
                        chainLabel.style.color = 'var(--primary-color)';
                        chainDiv.appendChild(chainLabel);
                        
                        const seqWrapper = document.createElement('div');
                        seqWrapper.style.display = 'flex';
                        seqWrapper.style.flexWrap = 'wrap';
                        seqWrapper.style.gap = '2px';
                        seqWrapper.style.fontFamily = 'monospace';
                        seqWrapper.style.fontSize = '1.1rem';
                        
                        nodesInChain.forEach(node => {
                            const dsspColor = dsspColors[node.dssp] || '#e0e0e0';
                            const aa = aminoAcidMap[node.residue] || '?';
                            const resSpan = document.createElement('span');
                            resSpan.textContent = aa;
                            resSpan.title = `${node.residue} ${node.position} (DSSP: ${node.dssp || 'Unknown'})`;
                            resSpan.style.cursor = 'pointer';
                            resSpan.style.padding = '2px 4px';
                            resSpan.style.borderRadius = '4px';
                            resSpan.style.transition = 'all 0.2s';
                            resSpan.style.backgroundColor = dsspColor;
                            resSpan.style.color = '#1f2937';
                            resSpan.style.fontWeight = 'bold';
                            
                            resSpan.onmouseover = () => {
                                resSpan.style.opacity = '0.7';
                            };
                            resSpan.onmouseout = () => {
                                resSpan.style.opacity = '1';
                            };
                            
                            resSpan.onclick = () => {
                                // Focus in 3D Force Graph
                                if(graph && node.x !== undefined) {
                                    const distance = 40;
                                    const distRatio = 1 + distance/Math.hypot(node.x, node.y, node.z);
                                    graph.cameraPosition(
                                        { x: node.x * distRatio, y: node.y * distRatio, z: node.z * distRatio },
                                        node,
                                        1000
                                    );
                                }
                                
                                // Focus in Molstar
                                if(molstarViewer) {
                                    molstarViewer.visual.clearSelection();
                                    molstarViewer.visual.select({ data: [{ struct_asym_id: node.chain, start_residue_number: parseInt(node.position), end_residue_number: parseInt(node.position) }] });
                                    molstarViewer.visual.focus([{ struct_asym_id: node.chain, start_residue_number: parseInt(node.position), end_residue_number: parseInt(node.position) }]);
                                }
                                
                                // Update interaction table
                                updateInteractionTable(node, graphData.links, metricsData);
                            };
                            
                            seqWrapper.appendChild(resSpan);
                        });
                        
                        chainDiv.appendChild(seqWrapper);
                        seqContent.appendChild(chainDiv);
                    });
                    
                    seqContainer.style.display = 'block';
                };
                
                renderSequence();
                
                // Inter-chain Interactions
                const renderInterchainTable = () => {
                    const container = document.getElementById('interchain-container');
                    const tbody = document.querySelector('#interchain-table tbody');
                    const noMsg = document.getElementById('no-interchain-msg');
                    const table = document.getElementById('interchain-table');
                    
                    tbody.innerHTML = '';
                    
                    const interchainLinks = graphData.links.filter(link => {
                        const s = typeof link.source === 'object' ? link.source : graphData.nodes.find(n => n.id === link.source);
                        const t = typeof link.target === 'object' ? link.target : graphData.nodes.find(n => n.id === link.target);
                        return s && t && s.chain && t.chain && s.chain !== t.chain;
                    });
                    
                    container.style.display = 'block';
                    
                    if (interchainLinks.length === 0) {
                        table.style.display = 'none';
                        noMsg.style.display = 'block';
                    } else {
                        table.style.display = 'table';
                        noMsg.style.display = 'none';
                        
                        interchainLinks.forEach(link => {
                            const sNode = typeof link.source === 'object' ? link.source : graphData.nodes.find(n => n.id === link.source);
                            const tNode = typeof link.target === 'object' ? link.target : graphData.nodes.find(n => n.id === link.target);
                            
                            const tr = document.createElement('tr');
                            tr.style.borderBottom = '1px solid var(--panel-border)';
                            tr.style.cursor = 'pointer';
                            tr.onmouseover = () => tr.style.backgroundColor = 'rgba(255,255,255,0.05)';
                            tr.onmouseout = () => tr.style.backgroundColor = 'transparent';
                            
                            tr.onclick = () => {
                                if (sNode) {
                                    const distance = 40;
                                    const distRatio = 1 + distance/Math.hypot(sNode.x || 0, sNode.y || 0, sNode.z || 0);
                                    if(graph && sNode.x !== undefined) {
                                        graph.cameraPosition(
                                            { x: sNode.x * distRatio, y: sNode.y * distRatio, z: sNode.z * distRatio },
                                            sNode,
                                            1000
                                        );
                                    }
                                    if(molstarViewer) {
                                        molstarViewer.visual.clearSelection();
                                        molstarViewer.visual.select({ data: [{ struct_asym_id: sNode.chain, start_residue_number: parseInt(sNode.position), end_residue_number: parseInt(sNode.position) }] });
                                        molstarViewer.visual.focus([{ struct_asym_id: sNode.chain, start_residue_number: parseInt(sNode.position), end_residue_number: parseInt(sNode.position) }]);
                                    }
                                    updateInteractionTable(sNode, graphData.links, metricsData);
                                }
                            };
                            
                            const dist = link.distance ? link.distance : 'N/A';
                            const ang = link.angle ? link.angle : 'N/A';
                            const interaction = link.interaction ? link.interaction : 'N/A';
                            
                            tr.innerHTML = `
                                <td style="padding: 8px; color: var(--primary-color);">${sNode.id}</td>
                                <td style="padding: 8px; color: var(--primary-color);">${tNode.id}</td>
                                <td style="padding: 8px;">${interaction}</td>
                                <td style="padding: 8px;">${dist}</td>
                                <td style="padding: 8px;">${ang}</td>
                            `;
                            tbody.appendChild(tr);
                        });
                    }
                };
                
                renderInterchainTable();
                
                // Ligand Interactions
                const renderLigandTable = () => {
                    const container = document.getElementById('ligand-container');
                    const tbody = document.querySelector('#ligand-table tbody');
                    const noMsg = document.getElementById('no-ligand-msg');
                    const table = document.getElementById('ligand-table');
                    
                    tbody.innerHTML = '';
                    
                    const isLigand = (resName) => {
                        const standardAA = {
                            'ALA':1, 'ARG':1, 'ASN':1, 'ASP':1, 'CYS':1,
                            'GLN':1, 'GLU':1, 'GLY':1, 'HIS':1, 'ILE':1,
                            'LEU':1, 'LYS':1, 'MET':1, 'PHE':1, 'PRO':1,
                            'SER':1, 'THR':1, 'TRP':1, 'TYR':1, 'VAL':1
                        };
                        return resName && resName.trim() !== '' && !standardAA[resName];
                    };
                    
                    const ligandLinks = graphData.links.filter(link => {
                        const s = typeof link.source === 'object' ? link.source : graphData.nodes.find(n => n.id === link.source);
                        const t = typeof link.target === 'object' ? link.target : graphData.nodes.find(n => n.id === link.target);
                        if (!s || !t) return false;
                        return isLigand(s.residue) || isLigand(t.residue);
                    });
                    
                    container.style.display = 'block';
                    
                    if (ligandLinks.length === 0) {
                        table.style.display = 'none';
                        noMsg.style.display = 'block';
                    } else {
                        table.style.display = 'table';
                        noMsg.style.display = 'none';
                        
                        ligandLinks.forEach(link => {
                            const sNode = typeof link.source === 'object' ? link.source : graphData.nodes.find(n => n.id === link.source);
                            const tNode = typeof link.target === 'object' ? link.target : graphData.nodes.find(n => n.id === link.target);
                            
                            const tr = document.createElement('tr');
                            tr.style.borderBottom = '1px solid var(--panel-border)';
                            tr.style.cursor = 'pointer';
                            tr.onmouseover = () => tr.style.backgroundColor = 'rgba(255,255,255,0.05)';
                            tr.onmouseout = () => tr.style.backgroundColor = 'transparent';
                            
                            const ligNode = isLigand(sNode.residue) ? sNode : tNode;
                            const targetNode = isLigand(sNode.residue) ? tNode : sNode;
                            
                            tr.onclick = () => {
                                if (ligNode) {
                                    const distance = 40;
                                    const distRatio = 1 + distance/Math.hypot(ligNode.x || 0, ligNode.y || 0, ligNode.z || 0);
                                    if(graph && ligNode.x !== undefined) {
                                        graph.cameraPosition(
                                            { x: ligNode.x * distRatio, y: ligNode.y * distRatio, z: ligNode.z * distRatio },
                                            ligNode,
                                            1000
                                        );
                                    }
                                    if(molstarViewer) {
                                        molstarViewer.visual.clearSelection();
                                        molstarViewer.visual.select({ data: [{ struct_asym_id: ligNode.chain, start_residue_number: parseInt(ligNode.position), end_residue_number: parseInt(ligNode.position) }] });
                                        molstarViewer.visual.focus([{ struct_asym_id: ligNode.chain, start_residue_number: parseInt(ligNode.position), end_residue_number: parseInt(ligNode.position) }]);
                                    }
                                    updateInteractionTable(ligNode, graphData.links, metricsData);
                                }
                            };
                            
                            const dist = link.distance ? link.distance : 'N/A';
                            const ang = link.angle ? link.angle : 'N/A';
                            const interaction = link.interaction ? link.interaction : 'N/A';
                            
                            tr.innerHTML = `
                                <td style="padding: 8px; color: var(--primary-color);">${ligNode.id}</td>
                                <td style="padding: 8px; color: var(--primary-color);">${targetNode.id}</td>
                                <td style="padding: 8px;">${interaction}</td>
                                <td style="padding: 8px;">${dist}</td>
                                <td style="padding: 8px;">${ang}</td>
                            `;
                            tbody.appendChild(tr);
                        });
                    }
                };
                
                renderLigandTable();
                
                // Complexity Analysis
                document.getElementById('complexity-analysis-container').style.display = 'block';
                renderComplexityAnalysis(graphData, metricsData);

            } catch (e) {
                console.error("Failed to load graph", e);
                graphViewerDiv.innerHTML = '<div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--error-color)">Failed to load graph</div>';
            }
        } else {
            graphViewerDiv.innerHTML = '<div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--text-muted)">Graph data not available</div>';
        }

        // Add resize listener
        window.addEventListener('resize', () => {
            if (graph) {
                graph.width(graphViewerDiv.clientWidth);
                graph.height(graphViewerDiv.clientHeight);
            }
        });
    }

    function parseNetworkData(nodesText, edgesText) {
        const nodes = [];
        const links = [];
        
        // Parse nodes
        const nodeLines = nodesText.trim().split('\n');
        if (nodeLines.length > 1) {
            const headers = nodeLines[0].split('\t');
            const idIdx = headers.indexOf('NodeId');
            const chainIdx = headers.indexOf('Chain');
            const posIdx = headers.indexOf('Position');
            const resIdx = headers.indexOf('Residue');
            const dsspIdx = headers.indexOf('Dssp');
            
            for (let i = 1; i < nodeLines.length; i++) {
                const cols = nodeLines[i].split('\t');
                if (cols.length > idIdx) {
                    nodes.push({
                        id: cols[idIdx],
                        chain: cols[chainIdx],
                        position: cols[posIdx],
                        residue: cols[resIdx],
                        dssp: dsspIdx > -1 ? cols[dsspIdx] : ''
                    });
                }
            }
        }
        
        // Parse edges
        const edgeLines = edgesText.trim().split('\n');
        if (edgeLines.length > 1) {
            const headers = edgeLines[0].split('\t');
            const id1Idx = headers.indexOf('NodeId1');
            const id2Idx = headers.indexOf('NodeId2');
            const interIdx = headers.indexOf('Interaction');
            const distIdx = headers.indexOf('Distance');
            const angleIdx = headers.indexOf('Angle');
            
            for (let i = 1; i < edgeLines.length; i++) {
                const cols = edgeLines[i].split('\t');
                if (cols.length > id2Idx) {
                    links.push({
                        source: cols[id1Idx],
                        target: cols[id2Idx],
                        interaction: cols[interIdx],
                        distance: distIdx > -1 ? cols[distIdx] : '',
                        angle: angleIdx > -1 ? cols[angleIdx] : ''
                    });
                }
            }
        }
        
        return { nodes, links };
    }

    function parseMetricsData(metricsText) {
        const metrics = {};
        if (!metricsText) return metrics;
        const lines = metricsText.trim().split('\n');
        if (lines.length > 1) {
            const headers = lines[0].split('\t');
            const idIdx = headers.indexOf('NodeId');
            const degIdx = headers.indexOf('Degree');
            const ccIdx = headers.indexOf('Clustering_Coefficient');
            const bcIdx = headers.indexOf('Betweenness_Centrality');
            const ecIdx = headers.indexOf('Eigenvector_Centrality');
            
            for (let i = 1; i < lines.length; i++) {
                const cols = lines[i].split('\t');
                if (cols.length > idIdx) {
                    metrics[cols[idIdx]] = {
                        degree: degIdx > -1 ? cols[degIdx] : 'N/A',
                        clustering: ccIdx > -1 && cols[ccIdx] !== 'N/A' ? parseFloat(cols[ccIdx]).toFixed(4) : 'N/A',
                        betweenness: bcIdx > -1 && cols[bcIdx] !== 'N/A' ? parseFloat(cols[bcIdx]).toFixed(4) : 'N/A',
                        eigenvector: ecIdx > -1 && cols[ecIdx] !== 'N/A' ? parseFloat(cols[ecIdx]).toFixed(4) : 'N/A'
                    };
                }
            }
        }
        return metrics;
    }

    function updateInteractionTable(selectedNode, links, metricsData) {
        const tableContainer = document.getElementById('interaction-details-container');
        const header = document.getElementById('interaction-header');
        const tbody = document.querySelector('#interaction-table tbody');
        
        tbody.innerHTML = '';
        
        const nodeLinks = links.filter(link => {
            const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
            const targetId = typeof link.target === 'object' ? link.target.id : link.target;
            return sourceId === selectedNode.id || targetId === selectedNode.id;
        });
        
        if (nodeLinks.length > 0) {
            document.getElementById('interaction-title').textContent = `Interactions for ${selectedNode.id}`;
            
            const metrics = metricsData[selectedNode.id];
            const metricsContainer = document.getElementById('node-metrics');
            if (metrics) {
                metricsContainer.innerHTML = `
                    <span><strong>Degree:</strong> ${metrics.degree}</span>
                    <span><strong>Clustering:</strong> ${metrics.clustering}</span>
                    <span><strong>Betweenness:</strong> ${metrics.betweenness}</span>
                    <span><strong>Eigenvector:</strong> ${metrics.eigenvector}</span>
                `;
            } else {
                metricsContainer.innerHTML = '';
            }
            
            tableContainer.style.display = 'block';
            
            nodeLinks.forEach(link => {
                const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
                const targetId = typeof link.target === 'object' ? link.target.id : link.target;
                
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--panel-border)';
                const otherNode = (typeof link.source === 'object' ? link.source.id : link.source) === selectedNode.id ? (typeof link.target === 'object' ? link.target.id : link.target) : (typeof link.source === 'object' ? link.source.id : link.source);
                const dist = link.distance ? link.distance : 'N/A';
                const ang = link.angle ? link.angle : 'N/A';
                const interaction = link.interaction ? link.interaction : 'N/A';
                
                tr.innerHTML = `
                    <td style="padding: 8px; color: var(--primary-color);">${otherNode}</td>
                    <td style="padding: 8px;">${interaction}</td>
                    <td style="padding: 8px;">${dist}</td>
                    <td style="padding: 8px;">${ang}</td>
                `;
                tbody.appendChild(tr);
            });
            
            // Scroll to the table
            tableContainer.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
        } else {
            tableContainer.style.display = 'none';
        }
    }

    function renderComplexityAnalysis(graphData, metricsData) {
        const nodes = graphData.nodes;
        const links = graphData.links;
        
        const numNodes = nodes.length;
        const numEdges = links.length;
        const density = numNodes > 1 ? (2 * numEdges) / (numNodes * (numNodes - 1)) : 0;
        
        let sumClustering = 0;
        let validClusteringCount = 0;
        for (const key in metricsData) {
            if (metricsData[key].clustering && metricsData[key].clustering !== 'N/A') {
                sumClustering += parseFloat(metricsData[key].clustering);
                validClusteringCount++;
            }
        }
        const avgClustering = validClusteringCount > 0 ? (sumClustering / validClusteringCount) : 0;
        
        const adj = {};
        nodes.forEach(n => adj[n.id] = []);
        links.forEach(l => {
            const u = typeof l.source === 'object' ? l.source.id : l.source;
            const v = typeof l.target === 'object' ? l.target.id : l.target;
            if(adj[u] && adj[v]) {
                adj[u].push(v);
                adj[v].push(u);
            }
        });
        
        const visited = new Set();
        const components = [];
        nodes.forEach(n => {
            if (!visited.has(n.id)) {
                const comp = [];
                const queue = [n.id];
                visited.add(n.id);
                let head = 0;
                while(head < queue.length) {
                    const u = queue[head++];
                    comp.push(u);
                    for(let v of adj[u]) {
                        if(!visited.has(v)) {
                            visited.add(v);
                            queue.push(v);
                        }
                    }
                }
                components.push(comp);
            }
        });
        
        let largestComp = [];
        for(let comp of components) {
            if(comp.length > largestComp.length) largestComp = comp;
        }
        
        let avgShortestPath = "N/A";
        let isLcc = components.length > 1;
        
        if (largestComp.length > 1) {
            let totalDist = 0;
            let pairs = 0;
            for (let i = 0; i < largestComp.length; i++) {
                const startNode = largestComp[i];
                const distances = {};
                const queue = [{node: startNode, dist: 0}];
                distances[startNode] = 0;
                let head = 0;
                while(head < queue.length) {
                    const current = queue[head++];
                    const u = current.node;
                    const d = current.dist;
                    for (let v of adj[u]) {
                        if (distances[v] === undefined) {
                            distances[v] = d + 1;
                            queue.push({node: v, dist: d + 1});
                            totalDist += d + 1;
                            pairs++;
                        }
                    }
                }
            }
            if (pairs > 0) {
                avgShortestPath = (totalDist / pairs).toFixed(4);
            }
        }
        
        document.getElementById('metric-nodes').textContent = numNodes;
        document.getElementById('metric-edges').textContent = numEdges;
        document.getElementById('metric-density').textContent = density.toFixed(4);
        document.getElementById('metric-clustering').textContent = avgClustering.toFixed(4);
        document.getElementById('metric-path').textContent = isLcc && avgShortestPath !== "N/A" ? `${avgShortestPath} (LCC)` : avgShortestPath;
        
        const degrees = nodes.map(n => n.degree || 1);
        const degreeCounts = {};
        degrees.forEach(d => { degreeCounts[d] = (degreeCounts[d] || 0) + 1; });
        
        const sortedDegrees = Object.keys(degreeCounts).map(Number).sort((a,b) => a - b);
        const counts = sortedDegrees.map(d => degreeCounts[d]);
        
        const isDarkMode = document.documentElement.getAttribute('data-theme') === 'dark';
        const fontColor = isDarkMode ? '#e5e7eb' : '#1f2937';
        const gridColor = isDarkMode ? '#374151' : '#e5e7eb';
        
        const renderPlots = () => {
            const isLogBar = document.getElementById('log-bar-chk').checked;
            const isLogScatter = document.getElementById('log-scatter-chk').checked;
            
            const layoutBase = {
                paper_bgcolor: 'rgba(0,0,0,0)',
                plot_bgcolor: 'rgba(0,0,0,0)',
                font: { color: fontColor, family: 'Inter, sans-serif' },
                margin: { t: 20, r: 20, l: 40, b: 40 }
            };
            
            const layoutBar = { ...layoutBase,
                xaxis: { title: 'Degree (k)', gridcolor: gridColor, type: isLogBar ? 'log' : 'linear' },
                yaxis: { title: 'Count P(k)', gridcolor: gridColor, type: isLogBar ? 'log' : 'linear' }
            };
            
            const layoutScatter = { ...layoutBase,
                xaxis: { title: 'Degree (k)', gridcolor: gridColor, type: isLogScatter ? 'log' : 'linear' },
                yaxis: { title: 'Count P(k)', gridcolor: gridColor, type: isLogScatter ? 'log' : 'linear' }
            };
            
            if (typeof Plotly !== 'undefined') {
                Plotly.newPlot('degree-bar-plot', [{
                    x: sortedDegrees,
                    y: counts,
                    type: 'bar',
                    marker: { color: '#6366f1' }
                }], layoutBar, {responsive: true});
                
                Plotly.newPlot('degree-scatter-plot', [{
                    x: sortedDegrees,
                    y: counts,
                    mode: 'markers+lines',
                    type: 'scatter',
                    marker: { color: '#10b981', size: 8 },
                    line: { color: '#10b981' }
                }], layoutScatter, {responsive: true});
            }
        };
        
        renderPlots();
        
        document.getElementById('log-bar-chk').onchange = renderPlots;
        document.getElementById('log-scatter-chk').onchange = renderPlots;
    }
});
