const fs = require('fs');

function parseNetworkData(nodesText, edgesText) {
    const nodes = [];
    const links = [];
    
    const nodeLines = nodesText.trim().split('\n');
    if (nodeLines.length > 1) {
        const headers = nodeLines[0].split('\t');
        const idIdx = headers.indexOf('NodeId');
        for (let i = 1; i < nodeLines.length; i++) {
            const cols = nodeLines[i].split('\t');
            if (cols.length > idIdx) {
                nodes.push({ id: cols[idIdx] });
            }
        }
    }
    
    const edgeLines = edgesText.trim().split('\n');
    if (edgeLines.length > 1) {
        const headers = edgeLines[0].split('\t');
        const id1Idx = headers.indexOf('NodeId1');
        const id2Idx = headers.indexOf('NodeId2');
        for (let i = 1; i < edgeLines.length; i++) {
            const cols = edgeLines[i].split('\t');
            if (cols.length > id2Idx) {
                links.push({
                    source: cols[id1Idx],
                    target: cols[id2Idx]
                });
            }
        }
    }
    
    return { nodes, links };
}

const nodesText = fs.readFileSync('4AKE.nodes', 'utf-8');
const edgesText = fs.readFileSync('4AKE.edges', 'utf-8');
const graphData = parseNetworkData(nodesText, edgesText);

const nodeIds = new Set(graphData.nodes.map(n => n.id));
let missing = 0;
for (const link of graphData.links) {
    if (!nodeIds.has(link.source)) {
        console.log("Missing source:", link.source);
        missing++;
    }
    if (!nodeIds.has(link.target)) {
        console.log("Missing target:", link.target);
        missing++;
    }
}
console.log("Total missing from edges vs nodes:", missing);
console.log("Total nodes:", graphData.nodes.length);
console.log("Total edges:", graphData.links.length);
