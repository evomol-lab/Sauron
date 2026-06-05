const fs = require('fs');

const nodesText = fs.readFileSync('1afw.nodes', 'utf-8');
const edgesText = fs.readFileSync('1afw.edges', 'utf-8');
const metricsText = fs.readFileSync('1afw_network_metrics.tsv', 'utf-8');

function parseNetworkData(nodesText, edgesText) {
    const nodes = [];
    const links = [];
    const nodeLines = nodesText.split('\n');
    const edgesLines = edgesText.split('\n');
    let hasHeader = nodeLines[0].toLowerCase().includes('nodeid');
    const startIdx = hasHeader ? 1 : 0;
    for (let i = startIdx; i < nodeLines.length; i++) {
        if (nodeLines[i].trim() === '') continue;
        const cols = nodeLines[i].split('\t');
        const id = cols[0];
        nodes.push({
            id: id,
            chain: cols[1] || '',
            position: cols[2] || '',
            residue: cols[3] || '',
            dssp: cols[4] || ''
        });
    }
    const eHeader = edgesLines[0].split('\t');
    const n1Idx = eHeader.findIndex(h => h.toLowerCase().includes('nodeid1'));
    const n2Idx = eHeader.findIndex(h => h.toLowerCase().includes('nodeid2'));
    const distIdx = eHeader.findIndex(h => h.toLowerCase().includes('distance'));
    const angIdx = eHeader.findIndex(h => h.toLowerCase().includes('angle'));
    const intIdx = eHeader.findIndex(h => h.toLowerCase().includes('interaction'));
    
    for (let i = 1; i < edgesLines.length; i++) {
        if (edgesLines[i].trim() === '') continue;
        const cols = edgesLines[i].split('\t');
        if (cols.length > Math.max(n1Idx, n2Idx)) {
            links.push({
                source: cols[n1Idx],
                target: cols[n2Idx],
                distance: distIdx > -1 ? cols[distIdx] : null,
                angle: angIdx > -1 ? cols[angIdx] : null,
                interaction: intIdx > -1 ? cols[intIdx].split(':')[0] : null
            });
        }
    }
    return { nodes, links };
}

function parseMetricsData(metricsText) {
    const metrics = {};
    const lines = metricsText.split('\n');
    if (lines.length > 0) {
        const header = lines[0].split('\t');
        const idIdx = header.findIndex(h => h.toLowerCase().includes('nodeid'));
        const degIdx = header.findIndex(h => h.toLowerCase().includes('degree'));
        const ccIdx = header.findIndex(h => h.toLowerCase().includes('clustering'));
        const bcIdx = header.findIndex(h => h.toLowerCase().includes('betweenness'));
        const ecIdx = header.findIndex(h => h.toLowerCase().includes('eigenvector'));
        
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

const graphData = parseNetworkData(nodesText, edgesText);
const metricsData = parseMetricsData(metricsText);

function renderComplexityAnalysis(graphData, metricsData) {
    const nodes = graphData.nodes;
    const links = graphData.links;
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
    console.log("Success! Avg Path:", avgShortestPath);
}

renderComplexityAnalysis(graphData, metricsData);

