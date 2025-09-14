import React from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Link } from "react-router-dom";
import "./Algorithms.css";

const algorithms = [
  {
    id: "scp",
    name: "Subgraph Coverage Pattern Mining (SCP)",
    num: 1,
    description:
      "Subgraph Coverage Pattern Mining (SCP) is a graph mining approach designed to uncover recurring structural motifs across multiple graphs. Instead of focusing only on individual frequent subgraphs, SCP emphasizes discovering patterns that cover a large fraction of the dataset, ensuring that the discovered motifs are representative of the global graph collection. This makes SCP especially suitable for analyzing biomedical networks, where subtle but widespread structural similarities can indicate meaningful biological mechanisms.",
    details:
      "SCP operates by systematically exploring subgraphs within collections of networks and identifying those that provide the highest overall coverage. A subgraph is considered a strong candidate if it recurs across multiple graphs in the dataset, thereby representing a common structural unit. By balancing coverage with overlap constraints, SCP prevents redundancy (detecting nearly identical subgraphs repeatedly) and ensures interpretability of results. In biomedical contexts, SCP allows researchers to uncover network motifs linked to disease progression, therapeutic mechanisms, and adverse drug interactions. Unlike traditional frequent subgraph mining, which may yield many redundant or overly specific results, SCP prioritizes patterns that summarize and represent the dataset as a whole.",
    researchApplications: [
  { 
    icon: "🔄", 
    title: "Drug Repurposing", 
    desc: "Identifying recurring drug interaction motifs across different diseases to suggest new therapeutic uses for existing drugs." 
  },
  { 
    icon: "⚠️", 
    title: "Adverse Reaction Clustering", 
    desc: "Detecting common substructures linked to adverse effects, helping to predict and prevent harmful drug interactions." 
  },
  { 
    icon: "🧬", 
    title: "Biomarker Discovery", 
    desc: "Recognizing motifs that frequently appear in disease-specific networks, aiding in the identification of diagnostic or prognostic biomarkers." 
  },
  { 
    icon: "🧪", 
    title: "Target Validation", 
    desc: "Validating candidate biological pathways or protein complexes by observing whether their substructures recur consistently across datasets." 
  },
  { 
    icon: "🧩", 
    title: "Comparative Network Analysis", 
    desc: "Comparing disease networks, drug combinations, or cell-type-specific networks to highlight shared biological mechanisms." 
  },
],

    keyParameters: [
      {
        name: "Minimum Support Threshold",
        description:
          "Specifies the minimum fraction of graphs in which a subgraph must appear to be considered significant. Controls how widely distributed a pattern must be.",
      },
      {
        name: "Maximum Allowed Overlap",
        description:
          "Limits redundancy by controlling how much two subgraphs are allowed to overlap. A lower value forces SCP to report more distinct patterns.",
      },
      {
        name: "Minimum Coverage Support",
        description:
          "Defines the minimum proportion of the overall dataset that must be covered by the discovered subgraphs. Ensures that reported motifs are globally representative.",
      },
    ],
    howItWorks: [
      {
        icon: "🧩",
        text: "Represent each dataset (e.g., drug interaction network, protein-protein network) as a labeled graph with nodes and edges.",
      },
      {
        icon: "🔍",
        text: "Systematically search for subgraphs that occur in multiple graphs, focusing on those with wide coverage across the dataset.",
      },
      {
        icon: "📊",
        text: "Apply thresholds and overlap constraints (Support, Coverage, Redundancy) to filter relevant subgraphs while avoiding duplicates.",
      },
      {
        icon: "🌐",
        text: "Highlight subgraphs that represent core structural motifs, such as recurring protein clusters or drug interaction modules.",
      },
    ],
  },
  {
    id: "tkg",
    name: "Top-k Frequent Subgraph Mining for Uncertain Graphs (TUSM)",
    num: 2,
    description:
      "Discovers frequent functional substructures in uncertain biological networks.",
    details:
      "Top-k Frequent Subgraph Mining for Uncertain Graphs focuses on finding the most frequent subgraphs in uncertain biological graphs. This allows detection of recurring biological motifs, protein complexes, and functional modules, crucial for understanding disease pathways.",
researchApplications: [
  {
    icon: "🧬",
    title: "Pathway Analysis",
    desc: "Examining molecular interaction patterns within biological pathways to uncover how drugs influence cellular processes and disease progression."
  },
  {
    icon: "🧩",
    title: "Disease Module Identification",
    desc: "Detecting tightly connected clusters of genes or proteins associated with specific diseases, aiding in biomarker discovery and precision medicine."
  },
  {
    icon: "💊",
    title: "Drug Mechanism Prediction",
    desc: "Inferring the underlying mechanisms of drug action by analyzing recurring network motifs, supporting rational drug design and repurposing."
  }
],

    keyParameters: [
      {
        name: "k",
        description:
          "Specifies how many of the most frequent subgraphs should be returned. Larger k increases coverage but also computational cost.",
      },
      {
        name: "Epsilon (ε)",
        description:
          "Controls error tolerance in probabilistic support estimation. Smaller ε yields more accurate results but requires more computation.",
      },
      {
        name: "Delta (δ)",
        description:
          "Confidence level for probabilistic support bounds. Lower δ means higher certainty at the cost of additional samples.",
      },
    ],
    howItWorks: [
      { icon: "📐", text: "Model biological network as an uncertain graph." },
      { icon: "⛏️", text: "Mine for the top-k most frequent subgraphs." },
      {
        icon: "📈",
        text: "Rank discovered patterns by their frequency and biological relevance.",
      },
    ],
  },
];


const Algorithms = () => {
  return (
    <div className="algorithms-container">
      <div className="container mx-auto px-4 py-8">
        {/* Header Section */}
        <div className="text-center mb-12">
          <h1 className="algorithms-title">Graph Mining Algorithms</h1>
          <p className="algorithms-subtitle">
            Explore algorithms for analyzing biological networks and drug discovery
          </p>
        </div>

        {/* Algorithm Tabs */}
        <Tabs defaultValue={algorithms[0].id} className="algorithms-tabs">
          <TabsList className="grid w-full grid-cols-2">
            {algorithms.map((algo) => (
              <TabsTrigger key={algo.id} value={algo.id} className="text-lg bg-purple-200 hover:bg-purple-500 hover:text-white">
                {algo.name}
              </TabsTrigger>
            ))}
          </TabsList>

          {algorithms.map((algo) => (
  <TabsContent key={algo.id} value={algo.id} className="mt-8">
    <Card className="algorithm-detail-card">
      <CardHeader>
        <CardTitle className="flex items-center justify-between">
          {algo.name}
          <Button
            variant="default"
            size="lg"
            className="bg-purple-600 hover:bg-purple-700 text-white text-lg"
            asChild
          >
            <Link to={`/tool${algo.num}`} className="text-lg">
              TRY IT OUT
            </Link>
          </Button>
        </CardTitle>
        <CardDescription className="text-lg ">
          {algo.description}
        </CardDescription>
      </CardHeader>

      <CardContent className="space-y-10">
        {/* Overview Section */}
        <div className="p-6 rounded-2xl bg-purple-50 border-l-4 border-purple-400 shadow-sm">
          <h4 className="detail-section-title mb-3">Overview</h4>
          <p className="text-foreground leading-relaxed">{algo.details}</p>
        </div>

        {/* Research Applications (Use Cases) */}
        <div>
          <h4 className="detail-section-title mb-5">Research Applications</h4>
          <div className="grid sm:grid-cols-2 lg:grid-cols-3 gap-6">
            {algo.researchApplications.map((item, idx) => (
  <div
    key={idx}
    className="p-5 rounded-xl shadow-md bg-gradient-to-br from-purple-700 to-purple-500 hover:shadow-lg transition flex flex-col items-center text-center"
  >
    <span className="text-2xl mb-2">{item.icon}</span>
    <h5 className="text-white font-semibold text-lg">{item.title}</h5>
    <p className="text-purple-100 mt-1">{item.desc}</p>
  </div>
))}
          </div>
        </div>

        {/* Key Parameters Section */}
        {algo.keyParameters && (
          <div>
            <h4 className="detail-section-title mb-5">Key Parameters</h4>
            <div className="grid sm:grid-cols-2 gap-6">
              {algo.keyParameters.map((param, idx) => (
                <div
                  key={idx}
                  className="p-5 rounded-xl bg-white border shadow-sm hover:shadow-md transition"
                >
                  <h5 className="font-semibold text-purple-700 mb-2 text-lg">
                    {param.name}
                  </h5>
                  <p className="text-foreground leading-relaxed">
                    {param.description}
                  </p>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* How It Works Section */}
        <div>
          <h4 className="detail-section-title mb-5">How It Works</h4>
          <ol className="relative border-l border-purple-300 space-y-6 ml-4">
            {algo.howItWorks.map((step, idx) => (
              <li key={idx} className="ml-4">
                <div className="absolute -left-2 w-4 h-4 bg-purple-500 rounded-full"></div>
                <div className="flex items-center gap-3">
                  <span className="text-xl">{step.icon}</span>
                  <p className="text-foreground leading-relaxed">{step.text}</p>
                </div>
              </li>
            ))}
          </ol>
        </div>
      </CardContent>
    </Card>
  </TabsContent>
))}

        </Tabs>
      </div>
    </div>
  );
};

export default Algorithms;
