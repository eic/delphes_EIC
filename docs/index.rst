.. Documentation master file

Delphes EIC
===========

.. toctree::
   :maxdepth: 2
   :caption: Getting started:

   install.md

   
.. toctree::
   :maxdepth: 2
   :caption: Links:

   GitHub <https://github.com/eic/delphes_EIC>
   Delphes <https://github.com/delphes/delphes>
   
   
.. Add main page text here

====
delphes_EIC
====

Welcome to the ``delphes_EIC`` project. This code aims to deliver a fast-simulation model of baseline or specific detectors for studies to support the Electron-Ion Collider (EIC) project. This code has been used in a few publications and papers and parts of it are available for citation through Zenodo:

- Charm jets as a probe for strangeness at the future Electron-Ion Collider. Miguel Arratia (UC Riverside amd Jefferson Lab), Yulia Furletova (Jefferson Lab), T. J. Hobbs (Southern Methodist U. and Jefferson Lab), Frederick Olness (Southern Methodist U.), Stephen J. Sekula (Southern Methodist U.). https://arxiv.org/abs/2006.12520 [hep-ph]
- A Delphes card for the EIC yellow-report detector. Miguel Arratia (UC Riverside amd Jefferson Lab) and Stephen J. Sekula (Southern Methodist U.). https://arxiv.org/abs/2103.06886 [physics.ins-det].  DOI: https://doi.org/10.5281/zenodo.4592887.
- Science Requirements and Detector Concepts for the Electron-Ion Collider: EIC Yellow Report.R. Abdul Khalek, A. Accardi, J. Adam, D. Adamiak, W. Akers et al. e-Print: https://arxiv.org/abs/2103.05419 [physics.ins-det]


.. image:: https://github.com/eic/delphes_EIC/blob/master/images/EICDetector_3D_CCDIS_CharmJet.png

.. image:: https://github.com/eic/delphes_EIC/raw/master/images/EICDetector_3D_CCDIS_CharmJet_DisplacedVtx.png

What's New?
----

* EIC PID code has been used to create IdentificationMaps for the mRICH, barrelDIRC, and dualRICH. No more external code needed!

Installation Instructions
----

#. Install PYTHIA8.3 following the instructions below, because code-patching is necessary to fix a bug in the DIS simulation.
#. Install Delphes3 following: https://github.com/delphes/delphes



EIC Yellow-Report Detector Models
----

The current model we recommend is ```delphes_card_allsilicon_3T.tcl```, which is based on detector studies from the EIC Yellow Report (https://arxiv.org/abs/2103.05419). This model incorporates all-silicon tracking, and EMCAL and HCAL. PID system responses are provided by efficiency maps based on EIC PID code (https://gitlab.com/preghenella/pid/). 

- Magnetic field: 1.5 T or 3.0 T
- Solenoid length: 2.0 m
- Tracker radius: 80 cm
- An HCAL and and EMCAL
- Particle ID systems: based on the mRICH, barrel DIRC, and dual RICH concepts articulated by the EIC community; implemented using IdentificationMaps

We currently simulte DIS using Pythia8 within Delphes. Again, detailed instructions for patching and installing it are below. The command file (ending in `.cmnd`) shown here is suitable for DIS at EIC. 

Generating Events
----

Run generation command:

.. code:: bash
   # PYTHIA8.3 and Delphes don't compile together, so we recommend a
   # 2-step process for now.
   # Generate HepMC output using PYTHIA8, then read into Delphes
   ./DelphesHepMC cards/delphes_card_allsilicon_3T.tcl out.root input.hepmc

You can see examples of analysis code in the Delphes page above

Visualizing Events 
----

Run the Delphes visualization script on your ROOT file, using the input Delphes card (TCL) file to help it structure the detector in the display:

.. code:: bash

   root -l examples/EventDisplay.C'("delphes_card_allsilicon_3T.tcl","out.root")'
 
The two examples shown here are for neutral-current and charged-current event for beam energies of 10 GeV electron on 100 GeV proton (63 GeV center-of-mass energy). 

EIC Collider Variations
----

Beam energy recommended benchmarking points are (the order is hadron on lepton):

* 275 on 18 GeV
* 100 on 10 GeV
* 100 on 5 GeV
* 41 on 5 GeV


The "SimpleAnalysis" framework
----

See the dedicated SimpleAnalysis_ Documentation.

.. _SimpleAnalysis: SimpleAnalysis/README.md
