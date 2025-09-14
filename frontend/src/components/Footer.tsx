import React from 'react';
import { Link } from 'react-router-dom';
import './Footer.css';

const Footer = () => {
  return (
    <footer className="biominer-footer">
      <div className="container mx-auto px-4 py-8">
        <div className="grid grid-cols-1 md:grid-cols-4 gap-8">
          <div className="footer-section">
            <h3 className="text-lg font-semibold text-primary mb-4">BioMiner</h3>
            <p className="text-muted-foreground text-sm">
              Advancing drug discovery through cutting-edge graph mining algorithms for medical researchers worldwide.
            </p>
          </div>
          
          <div className="footer-section">
            <h4 className="text-lg font-semibold text-foreground mb-3">Navigation</h4>
            <ul className="space-y-2">
              <li><Link to="/" className="footer-link">Home</Link></li>
              <li><Link to="/algorithms" className="footer-link">Algorithms</Link></li>
              <li><Link to="/about" className="footer-link">About Us</Link></li>
              <li><Link to="/contact" className="footer-link">Contact Us</Link></li>
            </ul>
          </div>
          
          <div className="footer-section">
            <h4 className="text-lg font-semibold text-foreground mb-3">Resources</h4>
            <ul className="space-y-2">
              <li><a href="#" className="footer-link">Documentation</a></li>
              <li><a href="#" className="footer-link">API Reference</a></li>
              <li><a href="#" className="footer-link">Tutorials</a></li>
              <li><a href="#" className="footer-link">Research Papers</a></li>
            </ul>
          </div>
          
          <div className="footer-section">
            <h4 className="text-lg font-semibold text-foreground mb-3">Support</h4>
            <ul className="space-y-2">
              <li><a href="#" className="footer-link">Help Center</a></li>
              <li><a href="#" className="footer-link">Community</a></li>
              <li><a href="#" className="footer-link">Bug Reports</a></li>
              <li><a href="#" className="footer-link">Feature Requests</a></li>
            </ul>
          </div>
        </div>
        
        <div className="border-t border-border mt-8 pt-6 text-center">
          <p className="text-sm text-muted-foreground">
            © 2025 BioMiner. All rights reserved. Built for advancing medical research.
          </p>
        </div>
      </div>
    </footer>
  );
};

export default Footer;