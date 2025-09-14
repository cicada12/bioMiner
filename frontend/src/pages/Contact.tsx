import React, { useState } from 'react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { useToast } from '@/hooks/use-toast';
import './Contact.css';

const Contact = () => {
  const { toast } = useToast();
  const [formData, setFormData] = useState({
    name: '',
    email: '',
    institution: '',
    inquiryType: '',
    subject: '',
    message: ''
  });

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    toast({
      title: "Message Sent!",
      description: "We'll get back to you within 24 hours.",
    });
    setFormData({
      name: '',
      email: '',
      institution: '',
      inquiryType: '',
      subject: '',
      message: ''
    });
  };

  const handleInputChange = (field: string, value: string) => {
    setFormData(prev => ({ ...prev, [field]: value }));
  };

  return (
    <div className="contact-container">
      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="contact-hero">
          <h1 className="contact-title">Get In Touch</h1>
          <p className="contact-subtitle">
            Have questions about BioMiner? Need support with your research? We're here to help.
          </p>
        </div>

        <div className="contact-content">
          {/* Contact Form */}
          <Card className="contact-form-card">
            <CardHeader>
              <CardTitle className='text-purple-700'>Send us a Message</CardTitle>
              <CardDescription>
                Fill out the form below and we'll get back to you as soon as possible.
              </CardDescription>
            </CardHeader>
            <CardContent>
              <form onSubmit={handleSubmit} className="contact-form">
                <div className="form-row">
                  <div className="form-field">
                    <Label htmlFor="name">Full Name *</Label>
                    <Input
                      id="name"
                      value={formData.name}
                      onChange={(e) => handleInputChange('name', e.target.value)}
                      required
                      placeholder="Dr. Jane Smith"
                    />
                  </div>
                  <div className="form-field">
                    <Label htmlFor="email">Email Address *</Label>
                    <Input
                      id="email"
                      type="email"
                      value={formData.email}
                      onChange={(e) => handleInputChange('email', e.target.value)}
                      required
                      placeholder="jane.smith@university.edu"
                    />
                  </div>
                </div>

                <div className="form-row">
                  <div className="form-field">
                    <Label htmlFor="institution">Institution/Organization</Label>
                    <Input
                      id="institution"
                      value={formData.institution}
                      onChange={(e) => handleInputChange('institution', e.target.value)}
                      placeholder="University of Research"
                    />
                  </div>
                  <div className="form-field">
                    <Label htmlFor="inquiryType">Inquiry Type *</Label>
                    <Select onValueChange={(value) => handleInputChange('inquiryType', value)}>
                      <SelectTrigger>
                        <SelectValue placeholder="Select inquiry type" />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="general">General Information</SelectItem>
                        <SelectItem value="technical">Technical Support</SelectItem>
                        <SelectItem value="collaboration">Research Collaboration</SelectItem>
                        <SelectItem value="licensing">Licensing & Partnerships</SelectItem>
                        <SelectItem value="bug">Bug Report</SelectItem>
                        <SelectItem value="feature">Feature Request</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>

                <div className="form-field">
                  <Label htmlFor="subject">Subject *</Label>
                  <Input
                    id="subject"
                    value={formData.subject}
                    onChange={(e) => handleInputChange('subject', e.target.value)}
                    required
                    placeholder="Brief description of your inquiry"
                  />
                </div>

                <div className="form-field">
                  <Label htmlFor="message">Message *</Label>
                  <Textarea
                    id="message"
                    value={formData.message}
                    onChange={(e) => handleInputChange('message', e.target.value)}
                    required
                    rows={6}
                    placeholder="Please provide details about your inquiry, research needs, or any specific questions you have about BioMiner..."
                  />
                </div>

                <Button type="submit" size="lg" className="submit-button text-white bg-purple-500 font-bold">
                  SEND MESSAGE
                </Button>
              </form>
            </CardContent>
          </Card>

          {/* Contact Information */}
          <div className="contact-info">
            <Card className="contact-info-card">
              <CardHeader>
                <CardTitle>Contact Information</CardTitle>
              </CardHeader>
              <CardContent className="space-y-6">
                <div className="contact-method">
                  <h4 className="contact-method-title">Email Support</h4>
                  <p className="contact-method-detail">support@biominer.org</p>
                  <p className="contact-method-description">
                    For technical support and general inquiries
                  </p>
                </div>

                <div className="contact-method">
                  <h4 className="contact-method-title">Research Partnerships</h4>
                  <p className="contact-method-detail">partnerships@biominer.org</p>
                  <p className="contact-method-description">
                    For collaboration opportunities and licensing
                  </p>
                </div>

                <div className="contact-method">
                  <h4 className="contact-method-title">Phone</h4>
                  <p className="contact-method-detail">+1 (555) 123-4567</p>
                  <p className="contact-method-description">
                    Monday - Friday, 9 AM - 5 PM EST
                  </p>
                </div>

                <div className="contact-method">
                  <h4 className="contact-method-title">Mailing Address</h4>
                  <p className="contact-method-detail">
                    BioMiner Research Lab<br />
                    XXXXX
                  </p>
                </div>
              </CardContent>
            </Card>

            
          </div>
        </div>
      </div>
    </div>
  );
};

export default Contact;