import pandas as pd
import requests
import plotly.graph_objects as go
import networkx as nx
from typing import Dict, List, Optional
import streamlit as st
import logging
import re
import py3Dmol
from Bio.Align import PairwiseAligner
import streamlit as st

# Sidebar navigation
st.sidebar.header("Navigation")
page = st.sidebar.radio("Select Page", ["Home Page", "Pathway Explorer", "Tutorial"])

if page == "Home Page":
    st.title("PathwaySync: Home Page")
   st.image("helix-6189400_1920.jpg", use_container_width=True, caption="Decoding Molecular Pathways")
    # About Section
    st.header("About Page")
    st.markdown("""
    Welcome to **PathwaySync**!

    This app helps you:
    - Explore molecular pathways from **KEGG** and **Reactome**
    - Visualize interactive 3D protein structures
    - Get AI-driven enzyme suggestions

    **Objectives:**
    - Provide an integrated platform for pathway exploration
    - Help researchers visualize and analyze molecular data
    - Enable easy data export in JSON or CSV

    **Features:**
    1) Interactive pathway graphs  
    2) Protein 3D structure viewer  
    3) AI enzyme suggestion tool  
    4) Data export options  
    5) Simple, user-friendly interface  
    """)

    st.markdown("---")

    # Team Section
    st.header("Team Page")
    st.image("./WhatsApp Image 2025-05-08 at 9.57.26 PM.jpeg", width=150, caption="Avanti Pandit")
    st.subheader("Avanti Pandit")
    st.write("""
    MSc Bioinformatics student passionate about exploring the interface between computational biology and molecular medicine.
    I developed **PathwaySync** to help researchers and students explore biological pathways with ease.
    """)

    st.markdown("""**Acknowledgement:**
    We would like to express our sincere gratitude to Dr. Kushagra Kashyap for his invaluable guidance, support, and encouragement throughout the development of this mini-project. His insights and mentorship were instrumental in shaping PathwaySync into a successful tool""")

elif page == "Tutorial":
    # Existing tutorial section
    st.header("PathwaySync Tutorial")
    st.markdown("""
    **Welcome to PathwaySync!**  
    This app integrates KEGG and Reactome pathways, visualizes protein structures, and suggests enzymes. Here's how to use it:
    1. **Pathway Explorer**: Select a database (KEGG or Reactome) and enter a pathway ID (e.g., 'hsa00010' for KEGG, 'R-HSA-1640170' for Reactome).
    2. View the interactive pathway graph; nodes are color-coded (blue for genes/events, green for compounds).
    3. Enter a PDB ID (e.g., '1A2C') to visualize 3D protein structures.
    4. Input a protein sequence for AI-driven enzyme suggestions.
    5. Use the sidebar to export results as JSON or CSV.
    **Tips**: Ensure valid IDs and check logs if errors occur.
    """)

else:
    # Existing Pathway Explorer code (keep your current code here)
    # pass  # Replace with your Pathway Explorer code block

    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # st.set_page_config(page_title="PathwaySync", layout="wide")

    @st.cache_data
    def fetch_kegg_pathway(pathway_id: str) -> Dict:
        """Fetch pathway data from KEGG API with retry and enhanced parsing."""
        pathway_id = pathway_id.strip().lower()
        if not re.match(r"^[a-z]{3}\d{5}$", pathway_id):
            error_msg = "Invalid KEGG pathway ID. Use format like 'hsa00010'."
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)
            return {"nodes": [], "edges": [], "metadata": {}}
        with requests.Session() as session:
            return fetch_kegg_data(pathway_id, session, logger, st, parse_kegg_response)

    def fetch_kegg_data(pathway_id, session, logger, st, parse_kegg_response):
        try:
            url = f"http://rest.kegg.jp/get/{pathway_id}"
            response = session.get(url, timeout=10)
            response.raise_for_status()
            return parse_kegg_response(response.text)
        except requests.RequestException as e:
            error_msg = f"Error fetching KEGG data for {pathway_id}: {str(e)}"
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)
            return {"nodes": [], "edges": [], "metadata": {}}

    def parse_kegg_response(text: str) -> Dict:
        current_section = ""
        meta = {}
        node_list = []
        edge_list = []

        for line in text.split("\n"):
            if line.startswith("ENTRY"):
                meta["id"] = line.split()[1]
                meta["type"] = line.split()[2] if len(line.split()) > 2 else "Pathway"
            elif line.startswith("NAME"):
                meta["name"] = line.replace("NAME", "").strip()
            elif line.startswith("ORGANISM"):
                meta["organism"] = line.replace("ORGANISM", "").strip()
            elif line.startswith("GENE"):
                if not current_section:
                    current_section = "GENE"
                gene_id = line.split()[1]
                node_list.append({"id": gene_id, "type": "gene", "label": gene_id})
            elif line.startswith("COMPOUND"):
                if not current_section:
                    current_section = "COMPOUND"
                compound_id = line.split()[1]
                node_list.append({"id": compound_id, "type": "compound", "label": compound_id})
            elif line.startswith("REACTION"):
                if not current_section:
                    current_section = "REACTION"
                parts = line.split()
                if len(parts) > 2:
                    sources = parts[1:]
                    for i in range(len(sources) - 1):
                        edge_list.append((sources[i], sources[i + 1]))
            else:
                if line.startswith(" ") and current_section:
                    item_id = line.split()[0]
                    node_list.append({"id": item_id, "type": current_section.lower(), "label": item_id})

            if line.strip() == "///":
                break

        for gene in [n for n in node_list if n["type"] == "gene"]:
            for compound in [n for n in node_list if n["type"] == "compound"]:
                edge_list.append((gene["id"], compound["id"]))

        return {"nodes": node_list, "edges": edge_list, "metadata": meta}

    @st.cache_data
    def fetch_reactome_pathway(pathway_id: str) -> Dict:
        if not re.match(r"^R-[A-Z]{3}-\d+(?:\.\d+)?$", pathway_id):
            error_msg = "Invalid Reactome pathway ID. Use format like 'R-HSA-1640170'."
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)
            return {"nodes": [], "edges": [], "metadata": {}}
        
        try:
            url = f"https://reactome.org/ContentService/data/pathway/{pathway_id}/containedEvents"
            headers = {"Accept": "application/json"}
            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()
            return parse_reactome_response(response.json(), pathway_id)
        except requests.RequestException as e:
            error_msg = f"Error fetching Reactome data for {pathway_id}: {str(e)}"
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)
            return {"nodes": [], "edges": [], "metadata": {}}

    def parse_reactome_response(data: List, pathway_id: str) -> Dict:
        nodes = []
        edges = []
        metadata = {"id": pathway_id, "name": "Unknown", "type": "Pathway"}

        for event in data:
            if "stId" in event and "displayName" in event:
                nodes.append({"id": event["stId"], "type": "event", "label": event["displayName"]})
                metadata["name"] = event.get("displayName", metadata["name"])
                if "related" in event:
                    for rel in event["related"]:
                        if "stId" in rel:
                            if (event["stId"], rel["stId"]) not in edges:
                                edges.append((event["stId"], rel["stId"]))
                genes = [n["id"] for n in nodes if n["type"] == "gene"]
                compounds = [n["id"] for n in nodes if n["type"] == "compound"]
                for gene in genes:
                    for compound in compounds:
                        edges.append((gene, compound))

        return {"nodes": nodes, "edges": edges, "metadata": metadata}

    def create_pathway_graph(data: Dict) -> go.Figure:
        G = nx.DiGraph()
        for node in data["nodes"]:
            G.add_node(node["id"], type=node["type"], label=node["label"])
        G.add_edges_from(data["edges"])

        pos = nx.spring_layout(G, seed=42)
        edge_x, edge_y = [], []
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y, line=dict(width=1, color="#888"), hoverinfo="none", mode="lines"
        )

        node_x, node_y, node_text, node_colors = [], [], [], []
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(G.nodes[node]["label"])
            node_type = G.nodes[node]["type"]
            color_map = {"event": "skyblue", "reaction": "orange", "compound": "lightgreen"}
            node_colors.append(color_map.get(node_type, "gray"))

        node_trace = go.Scatter(
            x=node_x, y=node_y, mode="markers+text", text=node_text, textposition="top center",
            marker=dict(size=12, color=node_colors), hoverinfo="text"
        )

        fig = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                showlegend=False, hovermode="closest", margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False), yaxis=dict(showgrid=False, zeroline=False),
                title=f"Pathway: {data['metadata']['name']} ({data['metadata']['id']})"
            )
        )
        return fig

    def fetch_pdb_structure(pdb_id: str) -> Optional[str]:
        if not re.match(r"^[0-9A-Za-z]{4}$", pdb_id):
            error_msg = "Invalid PDB ID. Use format like '1A2C'."
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)
            return None
        
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            error_msg = f"Error fetching PDB data for {pdb_id}: {str(e)}"
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)
            return None

    def visualize_3d_structure(pdb_data: str) -> None:
        try:
            view = py3Dmol.view(width=400, height=300)
            view.addModel(pdb_data, "pdb")
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo()
            st.components.v1.html(view._make_html(), width=400, height=300)
        except Exception as e:
            error_msg = f"Error visualizing 3D structure: {str(e)}"
            logger.error(error_msg)
            try:
                st.error(error_msg)
            except NameError:
                print(error_msg)

    def get_reference_enzymes():
        return {
            "TPP2_HUMAN": "MATATEEPFFHGLLPKKETGAASFLCRYPEYDGRGVILAVLDTGDVDPAGMQVTTDG",
            "TRYP1_HUMAN": "MKPTTLLALLGAAAAVSAAPPNDDKIVGGYTCAANSIPYQVSLNSGSHFCGGSLISEQWVVSAAHCYK",
            "CDK2_HUMAN": "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLD",
            "MY_ENZYME_1": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK",
            "MY_ENZYME_2": "MKADRLVFGGKDGGVGKSAKPGVVAGTVGGAGKTVTVLTDGGGRQILG"
        }

    def suggest_missing_enzyme(sequence):
        try:
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            reference_enzymes = get_reference_enzymes()
            lines = sequence.strip().splitlines()
            if lines and lines[0].startswith(">"):
                lines = lines[1:]
            clean_sequence = "".join(lines).replace(" ", "").strip()
            if not clean_sequence:
                return "No valid sequence found. Please check your input."

            best_match = None
            best_score = 0

            for name, ref_seq in reference_enzymes.items():
                score = aligner.score(clean_sequence, ref_seq)
                max_possible = max(len(clean_sequence), len(ref_seq))
                similarity = (score / max_possible) * 100 if max_possible > 0 else 0
                if similarity > best_score:
                    best_score = similarity
                    best_match = name

            if best_match:
                return f"Suggested enzyme: {best_match} with {best_score:.2f}% similarity"
            else:
                return "No similar enzyme found in the database."
        except Exception as e:
            error_msg = f"Error in sequence alignment: {str(e)}"
            logger.error(error_msg)
            return f"Error analyzing sequence: {str(e)}"

    # Streamlit UI
    st.title("PathwaySync: Unified Molecular Pathway Explorer")
    st.write("Explore molecular pathways from KEGG and Reactome, visualize 3D structures, and get AI-driven enzyme suggestions.")

    st.sidebar.header("Navigation")
    page = st.sidebar.radio("Select Page", ["Pathway Explorer", "Tutorial"])

    if page == "Tutorial":
        st.header("PathwaySync Tutorial")
        st.markdown("""
        **Welcome to PathwaySync!**  
        This app integrates KEGG and Reactome pathways, visualizes protein structures, and suggests enzymes. Here's how to use it:
        1. **Pathway Explorer**: Select a database (KEGG or Reactome) and enter a pathway ID (e.g., 'hsa00010' for KEGG, 'R-HSA-1640170' for Reactome).
        2. View the interactive pathway graph; nodes are color-coded (blue for genes/events, green for compounds).
        3. Enter a PDB ID (e.g., '1A2C') to visualize 3D protein structures.
        4. Input a protein sequence for AI-driven enzyme suggestions.
        5. Use the sidebar to export results as JSON or CSV.
        **Tips**: Ensure valid IDs and check logs if errors occur.
        """)
    else:
        st.header("Pathway Explorer")
        db_choice = st.selectbox("Select Database", ["KEGG", "Reactome"])
        placeholder = "hsa00010" if db_choice == "KEGG" else "R-HSA-1640170"
        pathway_id = st.text_input(f"Enter {db_choice} Pathway ID (e.g., {placeholder})", "")

        if st.button("Fetch Pathway"):
            if pathway_id:
                with st.spinner(f"Fetching {db_choice} pathway data..."):
                    pathway_data = (
                        fetch_kegg_pathway(pathway_id) if db_choice == "KEGG"
                        else fetch_reactome_pathway(pathway_id)
                    )
                    if pathway_data and "nodes" in pathway_data and pathway_data["nodes"]:
                        st.session_state["pathway_data"] = pathway_data
                        st.session_state["pathway_id"] = pathway_id
                        st.session_state["db_choice"] = db_choice
                        st.success(f"{db_choice} pathway retrieved: {pathway_data['metadata']['name']}")
                        st.plotly_chart(create_pathway_graph(pathway_data), use_container_width=True)
                    else:
                        st.session_state["pathway_data"] = None
                        st.error("No data found for the given pathway ID or fetch failed.")
            else:
                st.error("Please enter a valid pathway ID.")

        st.subheader("3D Structure Viewer")
        pdb_id = st.text_input("Enter PDB ID (e.g., 1A2C)", "")
        if st.button("Visualize Structure"):
            if pdb_id:
                with st.spinner("Fetching PDB structure..."):
                    pdb_data = fetch_pdb_structure(pdb_id)
                    if pdb_data:
                        st.success("Structure retrieved successfully!")
                        visualize_3d_structure(pdb_data)
                    else:
                        st.error("Unable to retrieve structure.")
            else:
                st.error("Please enter a valid PDB ID.")

        st.subheader("AI Enzyme Suggestion")
        sequence = st.text_area("Enter Protein Sequence (FASTA or plain text)", "")
        if st.button("Suggest Enzyme"):
            if sequence:
                with st.spinner("Analyzing sequence..."):
                    suggestion = suggest_missing_enzyme(sequence)
                    st.write(suggestion)
            else:
                st.error("Please enter a valid protein sequence.")

    # Sidebar export section
    st.sidebar.header("Export Options")
    export_format = st.sidebar.selectbox("Export Format", ["JSON", "CSV"])
    pathway_data = st.session_state.get("pathway_data")
    pathway_id = st.session_state.get("pathway_id", "pathway")
    db_choice = st.session_state.get("db_choice", "database")

    if st.sidebar.button("Export", key="export_button"):
        if not pathway_data:
            st.sidebar.warning("No pathway data available. Please fetch a pathway first.")
        elif "nodes" not in pathway_data or not pathway_data["nodes"]:
            st.sidebar.error("Pathway data is empty. Please check the pathway ID or your network connection.")
        else:
            df_nodes = pd.DataFrame(pathway_data["nodes"])
            if export_format == "JSON":
                st.sidebar.download_button(
                    label="Download JSON",
                    data=df_nodes.to_json(orient="records"),
                    file_name=f"{pathway_id}_{db_choice}_pathway.json",
                    mime="application/json"
                )
            elif export_format == "CSV":
                st.sidebar.download_button(
                    label="Download CSV",
                    data=df_nodes.to_csv(index=False),
                    file_name=f"{pathway_id}_{db_choice}_pathway.csv",
                    mime="text/csv"
                )

    st.markdown("---")
    st.markdown("Built with Streamlit by [Avanti Pandit].")
