import React from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import './About.css';

const citations = [
  {
    paper: "A model of graph transactional coverage patterns with applications to drug discovery",
    authors: "Reddy, A. S., Reddy, P. K., Mondal, A., & Priyakumar, U. D. (2021). HiPC, 21–30.",
    link: "https://doi.org/10.1109/hipc53243.2021.00016"
  },
  {
    paper: "Mining subgraph coverage patterns from graph transactions",
    authors: "Reddy, A. S., Reddy, P. K., Mondal, A., & Priyakumar, U. D. (2021b). Int. J. of Data Science and Analytics, 13(2), 105–121.",
    link: "https://doi.org/10.1007/s41060-021-00292-y"
  },
  {
    paper: "Frequent subgraph pattern mining on uncertain graph data",
    authors: "Zou, Z., Li, J., Gao, H., & Zhang, S. (2009). CIKM, 583–592.",
    link: "https://doi.org/10.1145/1645953.1646028"
  },
  {
    paper: "TKG: Efficient Mining of Top-K frequent subgraphs",
    authors: "Fournier-Viger, P., Cheng, C., Lin, J. C., Yun, U., & Kiran, R. U. (2019). Lecture Notes in Computer Science, 209–226.",
    link: "https://doi.org/10.1007/978-3-030-37188-3_13"
  }
];

const About = () => {
  return (
    <div className="about-container">
      <div className="container mx-auto px-4 py-8">
        {/* Hero Section */}
        <div className="about-hero">
          <h1 className="about-title">About BioMiner</h1>
          <p className="about-subtitle">
            Pioneering the future of drug discovery through graph mining technologies
          </p>
        </div>

        {/* Mission + Citations Section */}
        <section className="mission-section">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-start">
            
            {/* Left - Mission */}
<Card className="shadow-lg rounded-2xl bg-white p-6">
  <CardHeader>
    <CardTitle className="text-purple-700 text-2xl font-semibold">Our Mission</CardTitle>
  </CardHeader>
  <CardContent>
    <p className="mission-text text-foreground">
      At here, we are dedicated to accelerating medical breakthroughs by providing researchers 
      with cutting-edge graph mining tools. Our platform bridges the gap between complex 
      computational algorithms and practical drug discovery applications.
    </p>
    <p className="mission-text text-foreground mt-4">
      We believe that by making advanced network analysis accessible to medical researchers 
      worldwide, we can dramatically reduce the time and cost of bringing life-saving 
      treatments to patients.
    </p>
  </CardContent>
</Card>

            {/* Right - Citations Card */}
            <Card className="bg-gradient-to-br from-purple-700 to-purple-500 shadow-lg rounded-2xl text-white">
              <CardHeader>
                <CardTitle className="text-xl2 font-semibold">Key Citations</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                {citations.map((cite, idx) => (
                  <div key={idx} className="border-b border-purple-400/40 pb-3">
                    {/* Title */}
                    <h5 className="font-semibold">{cite.paper}</h5>
                    {/* Authors */}
                    <p className="text-purple-100 text-sm">{cite.authors}</p>
                    {/* DOI Link */}
                    <a
                      href={cite.link}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="text-purple-200 text-xs underline hover:text-white"
                    >
                      {cite.link}
                    </a>
                  </div>
                ))}
              </CardContent>
            </Card>

          </div>
        </section>
      </div>
    </div>
  );
};

export default About;
