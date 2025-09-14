import React, { useState, useRef } from 'react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Progress } from '@/components/ui/progress';
import JSZip from "jszip";
import { saveAs } from "file-saver";
import { Carousel } from 'react-responsive-carousel';
import 'react-responsive-carousel/lib/styles/carousel.min.css';
import './Tool_scp.css';

type Step = 'upload' | 'parameters' | 'results';

const Tool1 = () => {
  const [currentStep, setCurrentStep] = useState<Step>('upload');
  const [progress, setProgress] = useState(0);

  // file states
  const [uploadedFile, setUploadedFile] = useState<File | null>(null);
  const [fileText, setFileText] = useState('');
  const [fetchedText, setFetchedText] = useState(''); // for preset datasets
  const [dataset, setDataset] = useState('');

  const [error, setError] = useState('');

  // parameters
  const [minSupport, setMinSupport] = useState(0.5);
  const [maxOverlap, setMaxOverlap] = useState(0.5);
  const [minCoverage, setMinCoverage] = useState(0.5);

  // results + UI
  const [isRunning, setIsRunning] = useState(false);
  const [results, setResults] = useState<any>(null);
  const [images, setImages] = useState<string[]>([]);
  const [isCarouselVisible, setCarouselVisible] = useState(false);

  // file preview modal
  const [showFileModal, setShowFileModal] = useState(false);

  const dropRef = useRef<HTMLDivElement | null>(null);

  const backendURL = "http://127.0.0.1:8000";

  // ------------------------
  // FILE HANDLING
  // ------------------------
  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    if (file.type !== "text/plain" && !file.name.endsWith(".txt")) {
      setError("Only .txt files are accepted.");
      return;
    }

    setError('');
    setUploadedFile(file);
    setDataset(''); // clear preset selection
    setFetchedText('');
    const reader = new FileReader();
    reader.onload = (ev: ProgressEvent<FileReader>) => {
      setFileText(ev.target?.result as string || "");
    };
    reader.readAsText(file);
  };

  const handleDrop = (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    const file = e.dataTransfer.files?.[0];
    if (!file) return;

    if (file.type !== "text/plain" && !file.name.endsWith(".txt")) {
      setError("Only .txt files are accepted.");
      return;
    }

    setError('');
    setUploadedFile(file);
    setDataset('');
    setFetchedText('');
    const reader = new FileReader();
    reader.onload = (ev: ProgressEvent<FileReader>) => {
      setFileText(ev.target?.result as string || "");
    };
    reader.readAsText(file);
  };

  const handleViewFile = () => {
    if (!uploadedFile && !fetchedText) return;

    // if user selected a preset dataset, show fetchedText
    if (uploadedFile && uploadedFile.name.startsWith("graph_transactional_dataset") && fetchedText) {
      setFileText(fetchedText);
      setShowFileModal(true);
      return;
    }

    // otherwise if uploadedFile, read it for preview
    if (uploadedFile) {
      const reader = new FileReader();
      reader.onload = (ev: ProgressEvent<FileReader>) => {
        setFileText(ev.target?.result as string || "");
        setShowFileModal(true);
      };
      reader.onerror = () => {
        setError("Failed to read file for preview.");
      };
      reader.readAsText(uploadedFile);
    }
  };

  // upload endpoint (called when moving on from upload step)
  const handleSubmit = async () => {
    if (!uploadedFile) return;
    const formData = new FormData();
    formData.append("file_upload", uploadedFile);

    try {
      const response = await fetch(`${backendURL}/uploadfile/`, {
        method: "POST",
        body: formData,
      });
      if (!response.ok) throw new Error(`Upload failed (${response.status})`);
      const result = await response.json();
      console.log("Upload successful:", result);
    } catch (err) {
      console.error("Error uploading file:", err);
      // setError("File upload failed. See console for details.");
    }
  };

  // ------------------------
  // PARAMETERS endpoint
  // ------------------------
  const submitParameters = async () => {
    try {
      const response = await fetch(`${backendURL}/set-parameters`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          min_support: minSupport,
          max_overlap: maxOverlap,
          min_coverage: minCoverage,
        }),
      });
      if (!response.ok) throw new Error("Failed to submit parameters");
      const data = await response.json();
      console.log("Parameters submitted:", data);
    } catch (err) {
      console.error("Error submitting parameters:", err);
      setError("Failed to submit parameters.");
    }
  };

  // ------------------------
  // RUN + GET IMAGES
  // ------------------------
  const runAlgorithm = async () => {
    setIsRunning(true);
    try {
      const response = await fetch(`${backendURL}/run`, { method: "POST" });
      if (!response.ok) throw new Error("Algorithm run failed");
      const data = await response.json();
      console.log("Run response:", data);
      setResults(data);

      // fetch images (if any)
      const imgRes = await fetch(`${backendURL}/get-images`);
      const imgData = await imgRes.json();
      if (imgData.images) setImages(imgData.images);
    } catch (err) {
      console.error("Error running algorithm:", err);
      // setError("Run failed. See console for details.");
    } finally {
      setIsRunning(false);
      setCurrentStep('results');
      setProgress(100);
    }
  };

  // ------------------------
  // DOWNLOAD RESULTS (zip)
  // ------------------------
  const handleDownloadResults = async () => {
    if (!images.length) return;
    const zip = new JSZip();
    const folder = zip.folder("results");

    for (let i = 0; i < images.length; i++) {
      const url = `${backendURL}${images[i]}`; // images[i] already contains /images/...
      const response = await fetch(url);
      const blob = await response.blob();
      const arrayBuffer = await blob.arrayBuffer();
      folder?.file(`result_${i + 1}${blob.type === "image/png" ? ".png" : ".jpg"}`, arrayBuffer);
    }

    const content = await zip.generateAsync({ type: "blob" });
    saveAs(content, "results.zip");
  };

  // modal helpers
  const openCarouselModal = () => setCarouselVisible(true);
  const closeCarouselModal = () => setCarouselVisible(false);
  const closeFileModal = () => setShowFileModal(false);

  // ------------------------
  // STEP NAVIGATION
  // ------------------------
  const handleNext = async () => {
    const order: Step[] = ['upload', 'parameters', 'results'];
    const i = order.indexOf(currentStep);
    if (i < order.length - 1) {
      // perform actions for the current step before advancing
      if (currentStep === 'upload') {
        await handleSubmit();
      }
      if (currentStep === 'parameters') {
        await submitParameters();
      }

      const nextStep = order[i + 1];
      setCurrentStep(nextStep);
      setProgress((i + 1) * 33);
    }
  };

  const handlePrevious = () => {
    const order: Step[] = ['upload', 'parameters', 'results'];
    const i = order.indexOf(currentStep);
    if (i > 0) {
      setCurrentStep(order[i - 1]);
      setProgress((i - 1) * 33);
    }
  };

  return (
    <div className="tool-container">
      <div className="container mx-auto px-4 py-8">
        {/* header */}
        <div
  className="tool-header"
  style={{
    position: 'relative',
    padding: '3rem 2rem',
    textAlign: 'center',
    color: 'white',
    borderRadius: '1rem',
    overflow: 'hidden',
    backgroundImage: 'url(/pw.png)', // replace with your image path
    backgroundSize: 'cover',
    backgroundPosition: 'bottom',
  }}
>
  {/* Optional overlay for better readability */}
  <div
    style={{
      position: 'absolute',
      top: 0,
      left: 0,
      width: '100%',
      height: '100%',
      backgroundColor: 'rgba(69, 0, 96, 0.54)', // dark overlay
      zIndex: 0,
    }}
  ></div>

  {/* Title and subtitle */}
  <h1 style={{ position: 'relative', zIndex: 1, fontSize: '2.5rem', marginBottom: '0.5rem' }}>
    SCP Mining Analysis Tool
  </h1>
  <p style={{ position: 'relative', zIndex: 1, fontSize: '1.25rem' }}>
    Follow these steps to analyze your dataset
  </p>
</div>

        {/* progress bar with border */}
        <div className="steps-progress mb-6">
          <div className="border border-gray-400 rounded-lg p-2">
            <Progress value={progress} className="progress-bar" />
          </div>
        </div>

        {/* step cards */}
        <div className="steps-list grid grid-cols-1 md:grid-cols-3 gap-4 mb-8">
          {[
            { id: 'upload', title: 'Upload Dataset', description: 'Upload your .txt dataset' },
            { id: 'parameters', title: 'Set Parameters', description: 'Configure algorithm parameters' },
            { id: 'results', title: 'View Results', description: 'Check and review outputs' }
          ].map((s, idx) => (
            <Card
              key={s.id}
              className={`cursor-pointer transition-all ${currentStep === s.id ? 'ring-2 ring-purple-600 bg-purple-100' : 'hover:bg-gray-100'}`}
              onClick={() => {
                const order: Step[] = ['upload', 'parameters', 'results'];
                if (order.indexOf(s.id as Step) <= order.indexOf(currentStep)) {
                  setCurrentStep(s.id as Step);
                  setProgress(order.indexOf(s.id as Step) * 33);
                }
              }}
            >
              <CardHeader className="pb-2">
                <CardTitle className="text-base flex items-center gap-2">
                  <span className="h-6 w-6 flex items-center justify-center rounded-full bg-purple-600 text-white text-sm">{idx + 1}</span>
                  {s.title}
                </CardTitle>
                <CardDescription className="text-sm">{s.description}</CardDescription>
              </CardHeader>
            </Card>
          ))}
        </div>

        {/* step content */}
        <Card className="step-content-card">
          {/* UPLOAD */}
          {currentStep === 'upload' && (
            <div className="step-upload">
              <CardHeader>
                <CardTitle>Upload Dataset</CardTitle>
                <CardDescription>Only .txt files supported</CardDescription>
              </CardHeader>

              <CardContent className="drop-box space-y-6">
                {/* drag & drop area */}
                <div
                  className={`upload-box ${currentStep !== 'upload' ? 'disabled' : ''}`}
                  ref={dropRef}
                  onDrop={handleDrop}
                  onDragOver={(e) => e.preventDefault()}
                >
                  <label htmlFor="file-upload" className="file-upload-label cursor-pointer">
                    <div className="upload-icon text-3xl">☁️</div>
                    <div>
                      {uploadedFile ? (
                        <div className="file-display flex items-center justify-center gap-4">
                          <span className='text-black'>{uploadedFile.name}</span>
                          <Button variant="outline" className='text-white bg-purple-500' size="sm" onClick={handleViewFile} disabled={currentStep !== 'upload'}>
                            View Dataset
                          </Button>
                        </div>
                      ) : (
                        <span className="text-black">Drag and drop a txt file here, or click to select a txt file</span>
                      )}
                    </div>
                  </label>

                  {/* native file input (hidden) */}
                  <input
                    id="file-upload"
                    type="file"
                    accept=".txt"
                    onChange={handleFileChange}
                    className="file-input hidden"
                    disabled={currentStep !== 'upload'}
                  />

                  {/* dataset dropdown for presets */}
                  <div className="dropdown-group mt-4">
                    <span className="block mb-1 text-black">Or use available ones</span>
                    <select
                      value={dataset}
                      onChange={async (e) => {
                        const selected = e.target.value;
                        if (selected === '') {
                          setDataset('');
                          setUploadedFile(null);
                          setFetchedText('');
                          setFileText('');
                          return;
                        }

                        setDataset(selected);
                        try {
                          // NOTE: this assumes a relative /datasets/<file> available on frontend server
                          const response = await fetch(`/datasets/${selected}`);
                          const blob = await response.blob();
                          const file = new File([blob], selected, { type: blob.type || 'text/plain' });
                          setUploadedFile(file);
                          const text = await blob.text();
                          setFetchedText(text);
                          setFileText('');
                        } catch (err) {
                          console.error('Error loading dataset:', err);
                          setFetchedText('');
                          setFileText('');
                          setUploadedFile(null);
                          setError('Failed to load preset dataset.');
                        }
                      }}
                      className="dropdown border rounded px-2 py-1 text-black"
                      disabled={currentStep !== 'upload'}
                    >
                      <option value="">-- Select a dataset --</option>
                      <option value="graph_transactional_dataset1.txt">graph_transactional_dataset1.txt</option>
                      <option value="graph_transactional_dataset2.txt">graph_transactional_dataset2.txt</option>
                      <option value="graph_transactional_dataset3.txt">graph_transactional_dataset3.txt</option>
                      <option value="graph_transactional_dataset4.txt">graph_transactional_dataset4.txt</option>
                    </select>
                  </div>
                </div>

                {error && <p className="error-message text-red-500">{error}</p>}

                <Button onClick={handleNext} disabled={!uploadedFile} className="w-full bg-purple-600 text-white hover:bg-purple-700">
                  Continue to Parameters
                </Button>
              </CardContent>
            </div>
          )}

          {/* PARAMETERS */}
          {currentStep === 'parameters' && (
            <div className="step-parameters">
              <CardHeader>
                <CardTitle>Set Parameters</CardTitle>
                <CardDescription>
                  Configure the algorithm parameters below. Adjust the sliders or type values to set thresholds for the SCP analysis.
                </CardDescription>
              </CardHeader>

              <CardContent className="space-y-6">
                {/* Min Support */}
                <div className="parameter-wrapper space-y-2">
                  <Label>Minimum Support</Label>
                  <input
                    type="range"
                    min="0"
                    max="1"
                    step="0.01"
                    value={minSupport}
                    onChange={(e) => setMinSupport(Number(e.target.value))}
                    className="w-full"
                  />
                  <Input
                    type="number"
                    min="0"
                    max="1"
                    step="0.01"
                    value={minSupport}
                    onChange={(e) => setMinSupport(Number(e.target.value))}
                  />
                  <p className="slider-info text-sm text-gray-600">
                    Minimum percentage of graphs a subgraph must appear in to be considered frequent.
                  </p>
                </div>

      {/* Max Overlap */}
      <div className="parameter-wrapper space-y-2">
        <Label>Maximum Allowed Overlap</Label>
        <input
          type="range"
          min="0"
          max="1"
          step="0.01"
          value={maxOverlap}
          onChange={(e) => setMaxOverlap(Number(e.target.value))}
          className="w-full"
        />
        <Input
          type="number"
          min="0"
          max="1"
          step="0.01"
          value={maxOverlap}
          onChange={(e) => setMaxOverlap(Number(e.target.value))}
        />
        <p className="slider-info text-sm text-gray-600">
          Maximum percentage of overlap allowed between discovered subgraphs.
        </p>
      </div>

      {/* Min Coverage */}
      <div className="parameter-wrapper space-y-2">
        <Label>Minimum Coverage</Label>
        <input
          type="range"
          min="0"
          max="1"
          step="0.01"
          value={minCoverage}
          onChange={(e) => setMinCoverage(Number(e.target.value))}
          className="w-full"
        />
        <Input
          type="number"
          min="0"
          max="1"
          step="0.01"
          value={minCoverage}
          onChange={(e) => setMinCoverage(Number(e.target.value))}
        />
        <p className="slider-info text-sm text-gray-600">
          Minimum number of graphs that the selected subgraphs should collectively cover.
        </p>
      </div>

      <p className="guidelines text-sm text-gray-500 mt-4">
        Guidelines: Start with a high minimum support and minimum coverage, then gradually reduce them to explore more patterns. Begin with a low maximum overlap and increase it gradually if needed.
      </p>

      {/* Navigation */}
      <div className="navigation-buttons flex justify-between mt-6">
        <Button variant="outline" onClick={handlePrevious}>Previous</Button>
        <Button color="purple" onClick={async () => { 
          await submitParameters(); 
          await runAlgorithm(); 
        }}>
          Run Algorithm
        </Button>
      </div>
    </CardContent>
  </div>
)}

{/* RESULTS */}
{currentStep === 'results' && (
  <div className="step-results">
    <CardHeader>
      <CardTitle>Results</CardTitle>
    </CardHeader>
    <CardContent className="space-y-6">

      {/* Parameter Recap */}
      {minCoverage && (
        <div className="parameter-recap p-4 bg-purple-50 rounded-lg border-l-4 border-purple-400">
          <h4 className="font-semibold text-purple-700 mb-2">Analysis Parameters</h4>
          <p><strong>Min Coverage:</strong> {minCoverage}</p>
          <p><strong>Min Support:</strong> {minSupport}</p>
          <p><strong>Max Overlap:</strong> {maxOverlap}</p>
        </div>
      )}


      {/* Action Buttons */}
      <div className="flex flex-col sm:flex-row gap-4 mt-6">
        <Button onClick={openCarouselModal} className="bg-purple-600 text-white flex-1">
          View Results
        </Button>
        <Button onClick={handleDownloadResults} variant="outline" className="flex-1">
          Download Results
        </Button>
        <Button
          variant="outline"
          onClick={() => {
            setCurrentStep('upload');
            setProgress(0);
            setUploadedFile(null);
            setResults(null);
            setImages([]);
          }}
          className="flex-1"
        >
          Start New Analysis
        </Button>
      </div>
    </CardContent>

    {/* Carousel modal */}
    {isCarouselVisible && (
      <div className="modal-overlay" onClick={closeCarouselModal}>
        <div className="modal-result-content" onClick={(e) => e.stopPropagation()}>
          <button onClick={closeCarouselModal} className="close-modal-button">X</button>
          <p className="result-text">Subgraph Coverage Patterns</p>
          <p className="result-text">
            Min Coverage: {minCoverage} &nbsp; Min Support: {minSupport} &nbsp; Max Overlap: {maxOverlap}
          </p>
          <div className="step-box">
            <Carousel showThumbs={false} infiniteLoop useKeyboardArrows>
              {images.map((url, idx) => (
                <div key={idx}>
                  <img src={`${backendURL}${url}`} alt={`SCP ${idx}`} />
                </div>
              ))}
            </Carousel>
          </div>
        </div>
      </div>
    )}
  </div>
)}
        </Card>
      </div>

      {/* File preview modal */}
      {showFileModal && (
        <div className="modal-overlay" onClick={closeFileModal}>
          <div className="modal-result-content" onClick={(e) => e.stopPropagation()}>
            <div className="modal-header">
              <h3>Dataset Preview: {uploadedFile?.name || dataset}</h3>
              <button onClick={closeFileModal} className="close-button">X</button>
            </div>
            <pre className="modal-content" style={{ whiteSpace: 'pre-wrap', maxHeight: 400, overflow: 'auto' }}>
              {fileText}
            </pre>
          </div>
        </div>
      )}
    </div>
  );
};

export default Tool1;
