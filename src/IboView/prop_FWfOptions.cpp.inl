// this code is GENERATED. Do not change it directly!
void FWfOptions::InitProperties() {
   m_RunScf = true;
   m_RunIbba = true;
   m_Charge = 0;
   m_ExtraSpin = 0;
   m_ScfMethod = "Kohn-Sham (DFJX-RKS)";
   m_ScfOptions = "";
   m_Functional = "PBE";
   m_OrbBasis = "def2-TZVP";
   m_FitBasis = "univ-JFIT";
   m_LocMethod = "IBO (exponent 2)";
   m_OrbDivision = "As input wf";
   m_OrbDisplay = "Occupied only";
   m_AoType = "IAO (Sym Orth.)";
   m_BondOrderType = "IAO/ReNorm";
   m_LevelShiftClosed = -1.0;
   m_LevelShiftOpen = -0.5;
   m_ThrGrad = 1e-05;
   m_ThrEnergy = 1e-06;
   m_MaxIter = 4096;
   m_WorkSpaceMb = 200;
   m_NumThreads = 0;
}

char const *BoolToCstr(bool o);

QString FWfOptions::GetOptionsDesc(FOptionsDescType Type) const {
   FOptionsDesc desc(Type, "WfOptions", "doc.wf_options()");
   using std::abs;
   if (m_RunScf != true)
      desc.setp(-1, "run_scf", m_RunScf);
   if (m_RunIbba != true)
      desc.setp(-1, "run_ibba", m_RunIbba);
   if (m_Charge != 0)
      desc.setp(-1, "charge", m_Charge);
   if (m_ExtraSpin != 0)
      desc.setp(-1, "extra_spin", m_ExtraSpin);
   if (m_ScfMethod != "Kohn-Sham (DFJX-RKS)")
      desc.setp(-1, "scf_method", m_ScfMethod);
   if (m_ScfOptions != "")
      desc.setp(-1, "scf_options", m_ScfOptions);
   if (m_Functional != "PBE")
      desc.setp(-1, "functional", m_Functional);
   if (m_OrbBasis != "def2-TZVP")
      desc.setp(-1, "orb_basis", m_OrbBasis);
   if (m_FitBasis != "univ-JFIT")
      desc.setp(-1, "fit_basis", m_FitBasis);
   if (m_LocMethod != "IBO (exponent 2)")
      desc.setp(-1, "loc_method", m_LocMethod);
   if (m_OrbDivision != "As input wf")
      desc.setp(-1, "orb_division", m_OrbDivision);
   if (m_OrbDisplay != "Occupied only")
      desc.setp(-1, "orb_display", m_OrbDisplay);
   if (m_AoType != "IAO (Sym Orth.)")
      desc.setp(-1, "ao_type", m_AoType);
   if (m_BondOrderType != "IAO/ReNorm")
      desc.setp(-1, "bond_order_type", m_BondOrderType);
   if (abs(m_LevelShiftClosed - -1.0) > 1e-4)
      desc.setp(-1, "level_shift_closed", m_LevelShiftClosed);
   if (abs(m_LevelShiftOpen - -0.5) > 1e-4)
      desc.setp(-1, "level_shift_open", m_LevelShiftOpen);
   if (abs(m_ThrGrad - 1e-05) > 1e-4)
      desc.setp(-1, "thr_grad", m_ThrGrad);
   if (abs(m_ThrEnergy - 1e-06) > 1e-4)
      desc.setp(-1, "thr_energy", m_ThrEnergy);
   if (m_MaxIter != 4096)
      desc.setp(-1, "max_iter", m_MaxIter);
   if (m_WorkSpaceMb != 200)
      desc.setp(-1, "work_space_mb", m_WorkSpaceMb);
   if (m_NumThreads != 0)
      desc.setp(-1, "num_threads", m_NumThreads);
   return s2q(desc.str());
}

void FWfOptions::setp(QVariantMap const &vm) {
   for (QVariantMap::const_iterator it = vm.begin(); it != vm.end(); ++ it) {
      QByteArray propName = it.key().toUtf8();
      int iMetaProp = this->metaObject()->indexOfProperty(propName.constData());
      if (iMetaProp == -1) {
         IvNotify(NOTIFY_Warning, IvFmt("{0.ClassName}::setp(): attempted to set non-existent property '%1' (to '%2'). Ignored.", it.key(), it.value().toString()));
         continue;
      }
      // could use QObject::setProperty(name, variant), but that might create dynamic
      // properties, which we here do not want. And since we have the meta property index anyway...
      QMetaProperty const &metaProp = this->metaObject()->property(iMetaProp);
      if (!metaProp.write(this, it.value())) {
         IvNotify(NOTIFY_Warning, IvFmt("{0.ClassName}::setp(): failed to set static property '%1' (type: %2, index: %3) to '%4' (type: %5) via meta-property interface. Incompatible types?", it.key(), metaProp.typeName(), iMetaProp, it.value().toString(), it.value().typeName()));
         continue;
      }
   }
}

void FWfOptions::SetRunScf(bool o) {
   if (o != m_RunScf) {
      m_RunScf = o;
      emit RunScfChanged(m_RunScf);
   }
}

void FWfOptions::SetRunIbba(bool o) {
   if (o != m_RunIbba) {
      m_RunIbba = o;
      emit RunIbbaChanged(m_RunIbba);
   }
}

void FWfOptions::SetCharge(int o) {
   if (o < -32) o = -32;
   if (o > 32) o = 32;
   if (o != m_Charge) {
      m_Charge = o;
      emit ChargeChanged(m_Charge);
   }
}

void FWfOptions::SetExtraSpin(int o) {
   if (o < 0) o = 0;
   if (o > 32) o = 32;
   if (o != m_ExtraSpin) {
      m_ExtraSpin = o;
      emit ExtraSpinChanged(m_ExtraSpin);
   }
}

void FWfOptions::SetScfMethod(QString const &o) {
   if (o != m_ScfMethod) {
      m_ScfMethod = o;
      emit ScfMethodChanged(m_ScfMethod);
   }
}

void FWfOptions::SetScfOptions(QString const &o) {
   if (o != m_ScfOptions) {
      m_ScfOptions = o;
      emit ScfOptionsChanged(m_ScfOptions);
   }
}

void FWfOptions::SetFunctional(QString const &o) {
   if (o != m_Functional) {
      m_Functional = o;
      emit FunctionalChanged(m_Functional);
   }
}

void FWfOptions::SetOrbBasis(QString const &o) {
   if (o != m_OrbBasis) {
      m_OrbBasis = o;
      emit OrbBasisChanged(m_OrbBasis);
   }
}

void FWfOptions::SetFitBasis(QString const &o) {
   if (o != m_FitBasis) {
      m_FitBasis = o;
      emit FitBasisChanged(m_FitBasis);
   }
}

void FWfOptions::SetLocMethod(QString const &o) {
   if (o != m_LocMethod) {
      m_LocMethod = o;
      emit LocMethodChanged(m_LocMethod);
   }
}

void FWfOptions::SetOrbDivision(QString const &o) {
   if (o != m_OrbDivision) {
      m_OrbDivision = o;
      emit OrbDivisionChanged(m_OrbDivision);
   }
}

void FWfOptions::SetOrbDisplay(QString const &o) {
   if (o != m_OrbDisplay) {
      m_OrbDisplay = o;
      emit OrbDisplayChanged(m_OrbDisplay);
   }
}

void FWfOptions::SetAoType(QString const &o) {
   if (o != m_AoType) {
      m_AoType = o;
      emit AoTypeChanged(m_AoType);
   }
}

void FWfOptions::SetBondOrderType(QString const &o) {
   if (o != m_BondOrderType) {
      m_BondOrderType = o;
      emit BondOrderTypeChanged(m_BondOrderType);
   }
}

void FWfOptions::SetLevelShiftClosed(double o) {
   if (o < -32.0) o = -32.0;
   if (o > 32.0) o = 32.0;
   if (o != m_LevelShiftClosed) {
      m_LevelShiftClosed = o;
      emit LevelShiftClosedChanged(m_LevelShiftClosed);
   }
}

void FWfOptions::SetLevelShiftOpen(double o) {
   if (o < -32.0) o = -32.0;
   if (o > 32.0) o = 32.0;
   if (o != m_LevelShiftOpen) {
      m_LevelShiftOpen = o;
      emit LevelShiftOpenChanged(m_LevelShiftOpen);
   }
}

void FWfOptions::SetThrGrad(double o) {
   if (o < 1e-14) o = 1e-14;
   if (o > 100000000.0) o = 100000000.0;
   if (o != m_ThrGrad) {
      m_ThrGrad = o;
      emit ThrGradChanged(m_ThrGrad);
   }
}

void FWfOptions::SetThrEnergy(double o) {
   if (o < 1e-14) o = 1e-14;
   if (o > 100000000.0) o = 100000000.0;
   if (o != m_ThrEnergy) {
      m_ThrEnergy = o;
      emit ThrEnergyChanged(m_ThrEnergy);
   }
}

void FWfOptions::SetMaxIter(int o) {
   if (o < 0) o = 0;
   if (o != m_MaxIter) {
      m_MaxIter = o;
      emit MaxIterChanged(m_MaxIter);
   }
}

void FWfOptions::SetWorkSpaceMb(int o) {
   if (o < 0) o = 0;
   if (o > 128000) o = 128000;
   if (o != m_WorkSpaceMb) {
      m_WorkSpaceMb = o;
      emit WorkSpaceMbChanged(m_WorkSpaceMb);
   }
}

void FWfOptions::SetNumThreads(int o) {
   if (o < 0) o = 0;
   if (o > 128000) o = 128000;
   if (o != m_NumThreads) {
      m_NumThreads = o;
      emit NumThreadsChanged(m_NumThreads);
   }
}

void FWfOptions::CopyPropertiesFrom(FWfOptions const &other) {
   SetRunScf(other.GetRunScf());
   SetRunIbba(other.GetRunIbba());
   SetCharge(other.GetCharge());
   SetExtraSpin(other.GetExtraSpin());
   SetScfMethod(other.GetScfMethod());
   SetScfOptions(other.GetScfOptions());
   SetFunctional(other.GetFunctional());
   SetOrbBasis(other.GetOrbBasis());
   SetFitBasis(other.GetFitBasis());
   SetLocMethod(other.GetLocMethod());
   SetOrbDivision(other.GetOrbDivision());
   SetOrbDisplay(other.GetOrbDisplay());
   SetAoType(other.GetAoType());
   SetBondOrderType(other.GetBondOrderType());
   SetLevelShiftClosed(other.GetLevelShiftClosed());
   SetLevelShiftOpen(other.GetLevelShiftOpen());
   SetThrGrad(other.GetThrGrad());
   SetThrEnergy(other.GetThrEnergy());
   SetMaxIter(other.GetMaxIter());
   SetWorkSpaceMb(other.GetWorkSpaceMb());
   SetNumThreads(other.GetNumThreads());
}
// kate: syntax c++;
