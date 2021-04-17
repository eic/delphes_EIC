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

- <b><a href="https://inspirehep.net/literature/1802504">Charm jets as a probe for strangeness at the future Electron-Ion Collider</a></b>. <a href="https://inspirehep.net/authors/1203346">Miguel Arratia</a> (<a href="https://inspirehep.net/institutions/903304">UC, Riverside</a> and <a href="https://inspirehep.net/institutions/904961">Jefferson Lab</a>), Yulia Furletova (<a href="https://inspirehep.net/institutions/904961">Jefferson Lab</a>), <a href="https://inspirehep.net/authors/1057163">T.J. Hobbs</a> (<a href="https://inspirehep.net/institutions/905856">Southern Methodist U.</a> and <a href="https://inspirehep.net/institutions/904961">Jefferson Lab</a>), <a href="https://inspirehep.net/authors/994916">Fredrick Olness</a> (<a href="https://inspirehep.net/institutions/905856">Southern Methodist U.</a>), Stephen J. Sekula (<a href="https://inspirehep.net/institutions/905856">Southern Methodist U.</a>). e-Print: <a href="https://arxiv.org/abs/2006.12520">2006.12520</a>[hep-ph]
- <b><a href="https://inspirehep.net/literature/1851396">A Delphes card for the EIC yellow-report detector</a>.</b> <a href="https://inspirehep.net/authors/1203346">Miguel Arratia</a> (<a href="https://inspirehep.net/institutions/903304">UC, Riverside</a> and <a href="https://inspirehep.net/institutions/904961">Jefferson Lab</a>), <a href="https://inspirehep.net/authors/1021142">Stephen Sekula</a> (<a href="https://inspirehep.net/institutions/905856">Southern Methodist U.</a>). e-Print: <a href="https://arxiv.org/abs/2103.06886">2103.06886</a>[physics.ins-det]. DOI: <a href="https://doi.org/10.5281/zenodo.4592887">10.5281/zenodo.4592887</a>(publication)
- <b><a href="https://inspirehep.net/literature/1851258">Science Requirements and Detector Concepts for the Electron-Ion Collider: EIC Yellow Report</a></b>. <a href="https://inspirehep.net/authors/1706729">R. Abdul Khalek</a>, <a href="https://inspirehep.net/authors/1019006">A. Accardi</a>, J. Adam, D. Adamiak, W. Akers et al. e-Print: <a href="https://arxiv.org/abs/2103.05419"> 2103.05419 </a>[physics.ins-det]

<p>
<img style="float:left; width:400px;" src="https://github.com/eic/delphes_EIC/blob/master/images/EICDetector_3D_CCDIS_CharmJet.png"/>
<img style="float:right; width:400px;" src="https://github.com/eic/delphes_EIC/raw/master/images/EICDetector_3D_CCDIS_CharmJet_DisplacedVtx.png"/>
</p>

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

### Generating Events

Run generation command:

.. code:: bash

   ./DelphesPythia8 cards/delphes_card_allsilicon_3T.tcl CC_DIS.cmnd out.root

You can see examples of analysis code in the Delphes page above

### Visualizing Events 

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
