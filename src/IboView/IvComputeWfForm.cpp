/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#include "Iv.h"
#include "IvComputeWfForm.h"
#include "IvSettings.h"
#include "QPropertyModel.h"
#include "ui_ComputeWfForm.h"
#include "IvDocument.h"
#include "CtBasisSet.h"
#include "CtBasisDesc.h" // for iLikelyEcpCharge
#include "CtAtomSet.h"
#include "CtRhfOptions.h"

extern "C" {
   size_t getMemorySize(); // memory_size.c
}

FComputeWfForm::FComputeWfForm(FDocument *document, QWidget *parent)
   : QDialog(parent),
     ui(new Ui::ComputeWfForm),
     m_pDocument(document),
     m_bMemoryOkay(true),
     m_bOtherError(false)
{
   ui->setupUi(this);
//    ui->label_WorkSpace->setVisible(false);
//    ui->spinBox_WorkSpacePerThread->setVisible(false);
   ui->lineEdit_ScfOptions->setValidator(0);
   ui->spinBox_WorkSpacePerThread->setEnabled(false);

//    ui->tabWidget->layout()->setContentsMargins(0, 0, 0, 0); // left, top, right, bottom

//    m_pDocument->GetWfOptions()->setObjectName("WfOptions");
   LinkPropertyWidgets(m_pDocument->GetWfOptions(), this, "wf_option");
   connect(ui->checkBox_RunScf, SIGNAL(toggled(bool)), this, SLOT(ToggleScfPage(bool)));
   connect(ui->checkBox_RunIbba, SIGNAL(toggled(bool)), this, SLOT(ToggleIbbaPage(bool)));

   connect(ui->comboBox_ScfMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(RecomputeMemory()));
   connect(ui->comboBox_ScfOrbitalBasis, SIGNAL(currentIndexChanged(int)), this, SLOT(RecomputeMemory()));
   connect(ui->comboBox_ScfFitBasis, SIGNAL(currentIndexChanged(int)), this, SLOT(RecomputeMemory()));
   connect(ui->spinBox_NumThreads, SIGNAL(valueChanged(int)), this, SLOT(RecomputeMemory()));
   connect(ui->spinBox_WorkSpacePerThread, SIGNAL(valueChanged(int)), this, SLOT(RecomputeMemory()));
   connect(ui->spinBox_WfCharge, SIGNAL(valueChanged(int)), this, SLOT(RecomputeMemory()));
   connect(ui->spinBox_WfSpin, SIGNAL(valueChanged(int)), this, SLOT(RecomputeMemory()));
//    connect(ui->lineEdit_ScfOptions, SIGNAL(textChanged(QString)), this, SLOT(RecomputeMemory()));
   connect(ui->lineEdit_ScfOptions, SIGNAL(editingFinished(void)), this, SLOT(RecomputeMemory()));
   RecomputeMemory();

   CheckEcps();

   // if there already is electronic structure stuff, don't make new IBOs by default, unless
   // the user actually asked for it explicitly.
   if (m_pDocument->GetCurrentFrame() && m_pDocument->GetCurrentFrame()->HaveOrbitals()) {
      ui->checkBox_RunScf->setChecked(false);
   }

   ui->buttonBox->setFocus();

   if (!IvRestoreWindowSize("ComputeWindow/Size", this))
      IvGuessSubDialogSize(this);
}

// void FComputeWfForm::closeEvent(QCloseEvent *event)
// {
//    IvSaveWindowSize("ComputeWindow/Size", this);
//    return QDialog::closeEvent(event);
// }


FComputeWfForm::~FComputeWfForm()
{
   IvSaveWindowSize("ComputeWindow/Size", this);
   delete ui;
}

void FComputeWfForm::ToggleScfPage(bool Checked)
{
   ui->page_WfSetup->setEnabled(Checked);
   if (Checked)
      ui->tabWidget->setCurrentIndex(0);
   else
      ui->tabWidget->setCurrentIndex(1);
}

void FComputeWfForm::ToggleIbbaPage(bool Checked)
{
   ui->page_IbbaSetup->setEnabled(Checked);
   if (Checked)
      ui->tabWidget->setCurrentIndex(1);
}

void FComputeWfForm::SetMemoryText(QString s, FStatusClass Class)
{
   ui->label_MemoryGuess->setStyleSheet(QString("QLabel{padding: .2ex; font-size: 12pt; font-weight: bold; %1} QLabel::disabled{background: #444; color:#aaa}").arg(GetStatusStyle(Class)));
   ui->label_MemoryGuess->setText(s);
}


void FComputeWfForm::SetWfInfoText(QString s)
{
//    ui->label_MemoryGuess->setStyleSheet(QString("QLabel{padding: .2ex; font-size: 12pt; font-weight: bold; %1} QLabel::disabled{background: #444; color:#aaa}").arg(GetStatusStyle(Class)));
   ui->label_WfInfo->setText(s);
}


QString FormatSpinMultiplicity(int Ms2) {
   if (Ms2 < 0)
      return QString("invalid (ms2=%1)").arg(Ms2);
   switch(Ms2) {
      case 0: return "singlet";
      case 1: return "doublet";
      case 2: return "triplet";
      case 3: return "quartet";
      case 4: return "quintet";
      case 5: return "sextet";
      case 6: return "heptet";
      case 7: return "octet";
      default: return QString("%1-let").arg(Ms2+1);
   }
}


class FResourceEstimateError : public std::runtime_error
{
public:
   explicit FResourceEstimateError(std::string const &program, std::string const &msg)
      : std::runtime_error("resource estimate for " + program + ": " + msg), m_program(program), m_msg(msg)
   {}
   std::string const &msg() const { return m_msg; }
   std::string const &program() const { return m_program; }
protected:
   std::string
      m_program, m_msg;
};


struct FResourceRequirements {
   size_t
      // note: all given in *bytes*!
      nWorkSpaceShared,
      nMinWorkSpacePerThread,
      // memory *not* on the work space stacks
      nMemOnHeap; // memory *not* on the work space stacks
   FResourceRequirements();
};

FResourceRequirements::FResourceRequirements()
   : nWorkSpaceShared(0), nMinWorkSpacePerThread(0), nMemOnHeap(0)
{}


FResourceRequirements EstimateScfResourceRequirements(ct::FHfOptions const &ScfOpt, ct::FAtomSet const &Atoms, ct::FWfDecl const &WfDecl)
{
   // TODO:
   // - this stuff should probably be moved into MicroScf itself.
   //   Maybe into the SCF options, or some other sort of auxiliary object which
   //   actually knows the internals of the Fock builders.
   FResourceRequirements
      r;
   bool
      IsRestricted = ScfOpt.OrbType == ct::ORBTYPE_Restricted,
      IsOpenShell = !IsRestricted || WfDecl.nElecA() != WfDecl.nElecB();
   size_t
      // total number of occupied orbitals
      nOcc = 0;
   if (IsRestricted)
      nOcc = std::max(WfDecl.nElecA(), WfDecl.nElecB());
   else
      nOcc = WfDecl.nElecA() + WfDecl.nElecB();

   // instanciate the orbital and fitting sets to see how large they are.
   ct::FBasisSet
      OrbBasis(Atoms, ct::BASIS_Orbital);
   size_t
      nAo = OrbBasis.nFn();

   {
      size_t
         nMinWorkSpacePerThread_Jk = 0;
      bool
         // Does this requires an xc functional calculation? True for DFT methods.
         NeedXc = ScfOpt.UseXc(),
         // does it need an Hartree-Fock exchange matrix? True for HF itself or hybrid KS.
         NeedK = ScfOpt.UseExactExchange();
      bool
         // DF-J / DF-JX for pure functional DFT?
         DfjOnly = (ScfOpt.JkAlgo == ct::JKALGO_DfCoulombOnlyCache3ix ||
                  ScfOpt.JkAlgo == ct::JKALGO_DfCoulombOnly ||
                  (ScfOpt.JkAlgo == ct::JKALGO_DfAuto && !NeedK)),
         // DF-JK for Hartree-Fock or Hybrid DFT?
         Dfjk = (ScfOpt.JkAlgo == ct::JKALGO_DfCoulombAndExchange ||
               (ScfOpt.JkAlgo == ct::JKALGO_DfAuto && NeedK)),
         // integral-direct 4-index integrals? (no DF at all)
         Direct4ix = (ScfOpt.JkAlgo == ct::JKALGO_4ixDirect);

      assert(!(DfjOnly && Dfjk));
      if (DfjOnly) {
         ct::FBasisSet
            // instanciate Coulomb fitting (JFIT) basis.
            FitBasis(Atoms, ct::BASIS_JFit);
         size_t
            nFit = FitBasis.nFn();
         if (ScfOpt.JkAlgo == ct::JKALGO_DfCoulombOnlyCache3ix)
            // that is the estimate for the fully caching small-molecule 3ix DF-RKS.
            r.nMemOnHeap += sizeof(double) * (nFit*((nAo*(nAo+1))/2));
         nMinWorkSpacePerThread_Jk = sizeof(double)*(FitBasis.nFnOfLargestShell() * nAo*nAo);
         // ^- 2nd term: As intermediate for FormIntMNF while building the cached integrals.
         //    The latter could certainly be reduced if I find some time. But for the moment
         //    this will have to work...

         // density-fitting metric matrix (J) matrix decomposition
         r.nWorkSpaceShared += sizeof(double)*(nFit*nFit);
      } else if (Dfjk) {
         ct::FBasisSet
            // instanciate Coulomb+Exchange fitting (JKFIT) basis.
            FitBasis(Atoms, ct::BASIS_JkFit);
         size_t
            nFit = FitBasis.nFn();
         // density-fitting metric matrix (J) matrix decomposition
         r.nWorkSpaceShared += sizeof(double)*(nFit*nFit);
         // intermediate half-transformed (mu i|F) integrals for kc/ko in
         // MakeJk. Not finally decided on whether these should be on heap or
         // stack. Note: here nOcc = nClosed + nActive (restricted) or nOccA + nOccB (unrestriced).
         r.nMemOnHeap += sizeof(double)*(nAo * nFit * nOcc);
         // ^- 2nd term: intermediate for FormIntMNF in MakeJk
         nMinWorkSpacePerThread_Jk = sizeof(double)*(FitBasis.nFnOfLargestShell() * nAo*nAo);
      } else if (Direct4ix) {
         // this is our super-lame 4-ix integral direct Fock builder. It is
         // super slow, but does not need fitting sets and probably has the
         // lowest memory demands of exact-exchange capable Fock builders.
         // Also, it may be helpful for validation of intermediate results against
         // other programs. Not really meant for real-world use, though.

         // compressed density matrices and output j/k matrices in AccJk4i
         r.nWorkSpaceShared += sizeof(double)*(nAo*nAo*5);
         // accumulation buffers for triangular j and kc/ko matrices in AccJk4i
         nMinWorkSpacePerThread_Jk = (nAo+1)*(nAo+2)/2*(2 + IsOpenShell? 1 : 0);
         // (ab|cd) integral buffer for a single shell quartet inside AccJk4i
         nMinWorkSpacePerThread_Jk += sqr(sqr(OrbBasis.nFnOfLargestShell()));
      } else {
         IvNotify(NOTIFY_Warning, IvFmt("unrecognized JkAlgo %1 --- cannot compute required memory estimate", int(ScfOpt.JkAlgo)));
         throw FResourceEstimateError("SCF", fmt::format("JkAlgo {} resource estimate not programmed", int(ScfOpt.JkAlgo)));
      }

      // need xc integration space?
      size_t
         nWorkSpacePerThread_Dfti = 0;
      {
         if (NeedXc)
            nWorkSpacePerThread_Dfti = sizeof(double)*(2 * nAo*nAo);
            // ^-- 2 fock matrices per thread in AccXc
         else
            nWorkSpacePerThread_Dfti = 0;
      }
      r.nMinWorkSpacePerThread = std::max(nMinWorkSpacePerThread_Jk, nWorkSpacePerThread_Dfti);
      r.nMinWorkSpacePerThread += size_t(20)<<20;  // 20 MB
      r.nWorkSpaceShared += sizeof(double) * ((8 + (IsOpenShell?4:0))*nAo*nAo); // various Fock/density matrices. Didn't count.
   }

   return r;
}



void FComputeWfForm::RecomputeMemory()
{
   m_bMemoryOkay = true;
   m_bOtherError = true; // will be cleared at end of 'try{..}'
   try {
      FFrame
         *pFrame = m_pDocument->GetCurrentFrame();
      ct::FAtomSet
         *pOrigAtoms = 0;
      if (pFrame)
         pOrigAtoms = pFrame->pGetAtoms();
      if (pOrigAtoms == 0)
         return SetMemoryText("[no geometry]", STATUS_Confused);
      // make a copy of the atom set and assign it the currenty selected bases.
      ct::FAtomSet
         Atoms(*pOrigAtoms);
      ct::FHfOptions
         ScfOpt;
      ct::FWfDecl
         WfDecl;
      FWfOptions
         *pWfOptions = m_pDocument->GetWfOptions();
      pWfOptions->AssignBasisSets(&Atoms);
      pWfOptions->AssignScfOptions(ScfOpt);
      pWfOptions->AssignWfDecl(WfDecl, &Atoms);

      int
         nElecEcp = int(Atoms.nEcpElec()),
         nElecTotal = (Atoms.NuclearCharge() - pWfOptions->GetCharge());
      if (nElecTotal - nElecEcp < 0)
         return SetMemoryText("[nElec < 0]", STATUS_Confused);


      // report info for actual numbers of electrons and spin quantum number ms2
      {
         QString
            sElecEcp;
         if (nElecEcp != 0)
            sElecEcp = QString("%1 ecp, ").arg(nElecEcp);
         QString
            sWfInfo = QString("%1 electrons (%2%3%4, %5%6), %7 spin (ms2=%8)")
            .arg(nElecTotal + nElecEcp)
            .arg(sElecEcp)
            .arg(WfDecl.nElecA()).arg(QChar(0x03b1)) // <- unicode for greek letter alpha
            .arg(WfDecl.nElecB()).arg(QChar(0x03b2)) // <- unicode for greek letter beta
            .arg(FormatSpinMultiplicity(WfDecl.Ms2)).arg(WfDecl.Ms2);
         SetWfInfoText(sWfInfo);
      }

      size_t
         nThreads = (size_t)pWfOptions->GetNumThreads(),
         nSysMem = getMemorySize(); // amount of physical memory on the system.
      FResourceRequirements
         req = EstimateScfResourceRequirements(ScfOpt, Atoms, WfDecl);

      size_t
         // estimate for total amount of memory required
         nMemEstTotal = req.nMemOnHeap + req.nWorkSpaceShared + nThreads * req.nMinWorkSpacePerThread;

      pWfOptions->SetWorkSpaceMb((req.nMinWorkSpacePerThread + req.nWorkSpaceShared/nThreads)>>20);


      if (nSysMem == 0)
         nSysMem = size_t(8) << 30; // assume 8 GB.

//       IvEmit("This system has %1 GB of memory. Does it?", double(nSysMem)/double(1<<30));
      double
         f = (double(nMemEstTotal) + double(size_t(2)<<30)) / double(nSysMem);
      // convert to MB
      nMemEstTotal >>= 20;
      FStatusClass
         Class = STATUS_Confused;
      if (f < 0.3)
         Class = STATUS_Idle;
      else if (f < 0.8)
         Class = STATUS_Warning;
      else {
         Class = STATUS_Error;
         m_bMemoryOkay = false;
      }
      SetMemoryText(QString("%1 MB").arg((int)nMemEstTotal), Class);
      m_bOtherError = false;
   } catch (ct::FInputError const &e) {
      SetMemoryText("[input error]", STATUS_Confused);
      ui->textBrowser_WfNotes->setText(s2q(fmt::format("[Input error:]\n{}", e.what())));
      m_bOtherError = true;
   } catch (FResourceEstimateError const &e) {
      SetMemoryText(s2q("[" + e.msg() + "]"), STATUS_Confused);
      m_bOtherError = true;
   } catch (std::runtime_error &e) {
      SetMemoryText("[failed to load basis set]", STATUS_Confused);
      m_bOtherError = true;
   }
}

void FComputeWfForm::CheckEcps()
{
   try {
      FFrame
         *pFrame = m_pDocument->GetCurrentFrame();
      ct::FAtomSet
         *pOrigAtoms = 0;
      if (pFrame)
         pOrigAtoms = pFrame->pGetAtoms();
      if (pOrigAtoms == 0)
         return SetMemoryText("[no geometry]", STATUS_Confused);
      bool NeedsEcps = false;
      for (size_t iAt = 0; iAt < pOrigAtoms->size(); ++ iAt) {
         if ((*pOrigAtoms)[iAt].iElement > 36)
            NeedsEcps = true;
      }

      if (false && NeedsEcps) {
         ui->textBrowser_WfNotes->setText("<div style=\"font-size: 11pt; color:#ccc\"><p><div style=\"font-size: 13pt; color:white\">"
            "Sorry, MicroScf currently cannot do this calculation :(.</div> "
            "It can only compute Kohn-Sham wave functions up to element 36 (Kr), due to "
            "a lack of ECP integrals.</p>"
            "<p>Note: Chemical analysis of an imported wave function should still work.</p></div>");
         ui->checkBox_RunScf->setChecked(false);
         ui->checkBox_RunScf->setEnabled(false);
      }
   } catch (std::runtime_error const &e) {
      ui->textBrowser_WfNotes->setText("[failed to instanciate FAtomSet]");
   }
}





bool FComputeWfForm::GetRunScf() const
{
   return ui->checkBox_RunScf->isChecked();
}

bool FComputeWfForm::GetRunIbba() const
{
   return ui->checkBox_RunIbba->isChecked();
}

bool FComputeWfForm::IsMemoryOkay() const
{
   return m_bMemoryOkay;
}
