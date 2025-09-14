import React, { useState } from 'react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Link } from 'react-router-dom';
import './Home.css';
import { DropdownMenu, DropdownMenuContent, DropdownMenuItem, DropdownMenuTrigger } from "@/components/ui/dropdown-menu"

const algorithms = [
  {
    name: 'SCP Mining',
    id: 'Subgraph Coverage Pattern Mining',
    num: "1",
    description: 'Identifies influential nodes in protein-protein interaction networks',
    details: {
      researchApplications: [
        'Detect biologically significant motifs in drug or disease interaction networks',
        'Compare complex drug combinations through interaction structure analysis',
        'Identify shared mechanisms between treatments for different conditions',
      ],
    },
  },
  {
    name: 'TKG for Uncertain Graphs',
    id: 'Top k Frequent Subgraph Mining for Uncertain Graphs',
    num: "2",
    description: 'Discovers functional modules in biological networks',
    details: {
      researchApplications: [
        'Identify densely connected groups of proteins that function together',
        'Detect modules associated with specific diseases or biological pathways',
        'Explore how variations in network structure influence drug response',
      ],
    },
  },
];

const Home = () => {
  const [selectedAlgorithm, setSelectedAlgorithm] = useState(algorithms[0]);

  return (
    <div className="home-container">
      {/* Hero Section */}
      <section
        className="hero-section relative h-screen flex items-start justify-start bg-cover bg-top"
        style={{ backgroundImage: "url('../../public/hero1.png')" }}
      >
        {/* Overlay (optional) */}
        <div className="absolute inset-0 bg-black/40"></div>

        {/* Content */}
        <div className="relative z-10 text-left text-white px-8 pt-20 max-w-lg">
          <h1 className="hero-title text-4xl md:text-6xl font-bold">
            Accelerate Drug Discovery with
            <span className="separate"> Graph Mining</span>
          </h1>
          <p className="hero-subtitle mt-4 text-lg md:text-xl">
            Leverage cutting-edge graph algorithms to uncover hidden patterns in biological networks,
            identify novel drug targets, and advance medical research.
          </p>
          <div className="hero-actions mt-8 flex gap-4">
            <Button
              size="lg"
              asChild
              className="text-lg font-roboto bg-white text-purple-700 hover:shadow-lg hover:text-white transition-shadow"
            >
              <DropdownMenu>
                <DropdownMenuTrigger asChild>
                  <Button size="lg" variant="default" className="bg-purple-600 hover:bg-purple-700 text-white text-lg ">
                    Start Mining
                  </Button>
                </DropdownMenuTrigger>
                <DropdownMenuContent className="dropdown-menu-content bg-white shadow-lg rounded-xl p-2">
                  <DropdownMenuItem asChild className="dropdown-menu-item">
                    <Link to="/tool1" className="w-full px-3 py-2 rounded-md hover:bg-purple-100 font-semibold">
                      SCP
                    </Link>
                  </DropdownMenuItem>
                  <DropdownMenuItem asChild className="dropdown-menu-item">
                    <Link to="/tool2" className="w-full px-3 py-2 rounded-md hover:bg-purple-100 font-semibold">
                      TUSM
                    </Link>
                  </DropdownMenuItem>
                </DropdownMenuContent>
              </DropdownMenu>
            </Button>

            <Button
              size="lg"
              asChild
              className="text-lg font-roboto bg-white text-purple-700 hover:shadow-lg hover:text-white transition-shadow"
            >
              <Link to="/algorithms">Explore Algorithms</Link>
            </Button>
          </div>
        </div>
      </section>


{/* Algorithms Explorer Section */}
<section className="algorithms-explorer">
  <div className="container mx-auto px-4 py-16">
    <div className="text-left mb-12">
      <h2 className="section-title">Explore Our Algorithms</h2>
    </div>

    {/* Two-pane layout */}
<div className="grid grid-cols-1 md:grid-cols-4 gap-8">
  {/* Left Pane - Algorithm List */}
  <div className="algorithms-list md:col-span-1 space-y-4 border-r border-gray-300 pr-4">
    {algorithms.map((algo) => (
      <Card
        key={algo.id}
        className={`cursor-pointer transition transform hover:scale-[1.02] hover:shadow-md ${
          selectedAlgorithm.id === algo.id
            ? "ring-2 ring-purple-500 bg-purple-100"
            : "bg-white"
        }`}
        onClick={() => setSelectedAlgorithm(algo)}
      >
        <CardHeader className="pb-2">
          <CardTitle className="text-md font-semibold">{algo.name}</CardTitle>
          <CardDescription className="text-sm text-foreground line-clamp-2">
            {algo.description}
          </CardDescription>
        </CardHeader>
      </Card>
    ))}
  </div>

  {/* Right Pane - Algorithm Details */}
  <div className="algorithm-details md:col-span-3 pl-4">
    <Card className="h-full bg-white shadow-md rounded-xl">
      <CardHeader>
        <CardTitle className="flex items-center justify-between text-xl font-bold">
          {selectedAlgorithm.name}
          <Button className="bg-purple-600 text-white hover:bg-purple-700">
            <Link to={`/tool${selectedAlgorithm.num}`}>Try It Out</Link>
          </Button>
        </CardTitle>
      </CardHeader>

      <CardContent className="space-y-8">
        {/* Overview */}
        <div className="p-4 rounded-lg bg-purple-50 border-l-4 border-purple-400">
          <h4 className="font-semibold mb-2 text-purple-800 text-lg">Overview</h4>
          <p className="text-foreground">{selectedAlgorithm.description}</p>
        </div>

        {/* Research Applications */}
        <div>
          <h4 className="font-semibold text-lg mb-4">Research Applications</h4>
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            {selectedAlgorithm.details.researchApplications.map((item, idx) => (
              <div
                key={idx}
                className="p-3 rounded-lg bg-purple-500 text-white shadow-sm hover:shadow-md transition"
              >
                {item}
              </div>
            ))}
          </div>
        </div>
      </CardContent>
    </Card>
  </div>
</div>

  </div>
</section>


    </div>
  );
};

export default Home;