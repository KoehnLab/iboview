// this code is GENERATED. Do not change it directly!
   // update properties in *this from entries in variant map
   Q_SLOT void setp(QVariantMap const &vm);

   Q_PROPERTY(bool run_scf READ GetRunScf WRITE SetRunScf NOTIFY RunScfChanged)
   bool m_RunScf;
   Q_SLOT void SetRunScf(bool o);
   Q_SIGNAL void RunScfChanged(bool const &o);
   bool GetRunScf() const { return m_RunScf; };

   Q_PROPERTY(bool run_ibba READ GetRunIbba WRITE SetRunIbba NOTIFY RunIbbaChanged)
   bool m_RunIbba;
   Q_SLOT void SetRunIbba(bool o);
   Q_SIGNAL void RunIbbaChanged(bool const &o);
   bool GetRunIbba() const { return m_RunIbba; };

   Q_PROPERTY(int charge READ GetCharge WRITE SetCharge NOTIFY ChargeChanged)
   int m_Charge;
   Q_SLOT void SetCharge(int o);
   Q_SIGNAL void ChargeChanged(int const &o);
   int GetCharge() const { return m_Charge; };

   Q_PROPERTY(int extra_spin READ GetExtraSpin WRITE SetExtraSpin NOTIFY ExtraSpinChanged)
   int m_ExtraSpin;
   Q_SLOT void SetExtraSpin(int o);
   Q_SIGNAL void ExtraSpinChanged(int const &o);
   int GetExtraSpin() const { return m_ExtraSpin; };

   Q_PROPERTY(QString scf_method READ GetScfMethod WRITE SetScfMethod NOTIFY ScfMethodChanged)
   QString m_ScfMethod;
   Q_SLOT void SetScfMethod(QString const &o);
   Q_SIGNAL void ScfMethodChanged(QString const &o);
   QString const &GetScfMethod() const { return m_ScfMethod; };

   Q_PROPERTY(QString scf_options READ GetScfOptions WRITE SetScfOptions NOTIFY ScfOptionsChanged)
   QString m_ScfOptions;
   Q_SLOT void SetScfOptions(QString const &o);
   Q_SIGNAL void ScfOptionsChanged(QString const &o);
   QString const &GetScfOptions() const { return m_ScfOptions; };

   Q_PROPERTY(QString functional READ GetFunctional WRITE SetFunctional NOTIFY FunctionalChanged)
   QString m_Functional;
   Q_SLOT void SetFunctional(QString const &o);
   Q_SIGNAL void FunctionalChanged(QString const &o);
   QString const &GetFunctional() const { return m_Functional; };

   Q_PROPERTY(QString orb_basis READ GetOrbBasis WRITE SetOrbBasis NOTIFY OrbBasisChanged)
   QString m_OrbBasis;
   Q_SLOT void SetOrbBasis(QString const &o);
   Q_SIGNAL void OrbBasisChanged(QString const &o);
   QString const &GetOrbBasis() const { return m_OrbBasis; };

   Q_PROPERTY(QString fit_basis READ GetFitBasis WRITE SetFitBasis NOTIFY FitBasisChanged)
   QString m_FitBasis;
   Q_SLOT void SetFitBasis(QString const &o);
   Q_SIGNAL void FitBasisChanged(QString const &o);
   QString const &GetFitBasis() const { return m_FitBasis; };

   Q_PROPERTY(QString loc_method READ GetLocMethod WRITE SetLocMethod NOTIFY LocMethodChanged)
   QString m_LocMethod;
   Q_SLOT void SetLocMethod(QString const &o);
   Q_SIGNAL void LocMethodChanged(QString const &o);
   QString const &GetLocMethod() const { return m_LocMethod; };

   Q_PROPERTY(QString orb_division READ GetOrbDivision WRITE SetOrbDivision NOTIFY OrbDivisionChanged)
   QString m_OrbDivision;
   Q_SLOT void SetOrbDivision(QString const &o);
   Q_SIGNAL void OrbDivisionChanged(QString const &o);
   QString const &GetOrbDivision() const { return m_OrbDivision; };

   Q_PROPERTY(QString orb_display READ GetOrbDisplay WRITE SetOrbDisplay NOTIFY OrbDisplayChanged)
   QString m_OrbDisplay;
   Q_SLOT void SetOrbDisplay(QString const &o);
   Q_SIGNAL void OrbDisplayChanged(QString const &o);
   QString const &GetOrbDisplay() const { return m_OrbDisplay; };

   Q_PROPERTY(QString ao_type READ GetAoType WRITE SetAoType NOTIFY AoTypeChanged)
   QString m_AoType;
   Q_SLOT void SetAoType(QString const &o);
   Q_SIGNAL void AoTypeChanged(QString const &o);
   QString const &GetAoType() const { return m_AoType; };

   Q_PROPERTY(QString bond_order_type READ GetBondOrderType WRITE SetBondOrderType NOTIFY BondOrderTypeChanged)
   QString m_BondOrderType;
   Q_SLOT void SetBondOrderType(QString const &o);
   Q_SIGNAL void BondOrderTypeChanged(QString const &o);
   QString const &GetBondOrderType() const { return m_BondOrderType; };

   Q_PROPERTY(double level_shift_closed READ GetLevelShiftClosed WRITE SetLevelShiftClosed NOTIFY LevelShiftClosedChanged)
   double m_LevelShiftClosed;
   Q_SLOT void SetLevelShiftClosed(double o);
   Q_SIGNAL void LevelShiftClosedChanged(double const &o);
   double GetLevelShiftClosed() const { return m_LevelShiftClosed; };

   Q_PROPERTY(double level_shift_open READ GetLevelShiftOpen WRITE SetLevelShiftOpen NOTIFY LevelShiftOpenChanged)
   double m_LevelShiftOpen;
   Q_SLOT void SetLevelShiftOpen(double o);
   Q_SIGNAL void LevelShiftOpenChanged(double const &o);
   double GetLevelShiftOpen() const { return m_LevelShiftOpen; };

   Q_PROPERTY(double thr_grad READ GetThrGrad WRITE SetThrGrad NOTIFY ThrGradChanged)
   double m_ThrGrad;
   Q_SLOT void SetThrGrad(double o);
   Q_SIGNAL void ThrGradChanged(double const &o);
   double GetThrGrad() const { return m_ThrGrad; };

   Q_PROPERTY(double thr_energy READ GetThrEnergy WRITE SetThrEnergy NOTIFY ThrEnergyChanged)
   double m_ThrEnergy;
   Q_SLOT void SetThrEnergy(double o);
   Q_SIGNAL void ThrEnergyChanged(double const &o);
   double GetThrEnergy() const { return m_ThrEnergy; };

   Q_PROPERTY(int max_iter READ GetMaxIter WRITE SetMaxIter NOTIFY MaxIterChanged)
   int m_MaxIter;
   Q_SLOT void SetMaxIter(int o);
   Q_SIGNAL void MaxIterChanged(int const &o);
   int GetMaxIter() const { return m_MaxIter; };

   Q_PROPERTY(int work_space_mb READ GetWorkSpaceMb WRITE SetWorkSpaceMb NOTIFY WorkSpaceMbChanged)
   int m_WorkSpaceMb;
   Q_SLOT void SetWorkSpaceMb(int o);
   Q_SIGNAL void WorkSpaceMbChanged(int const &o);
   int GetWorkSpaceMb() const { return m_WorkSpaceMb; };

   Q_PROPERTY(int num_threads READ GetNumThreads WRITE SetNumThreads NOTIFY NumThreadsChanged)
   int m_NumThreads;
   Q_SLOT void SetNumThreads(int o);
   Q_SIGNAL void NumThreadsChanged(int const &o);
   int GetNumThreads() const { return m_NumThreads; };

   void InitProperties(); // set default values of properties and emit corresponding signals.
   // make script code for setting/collecting options in *this which were changed from their defaults
   QString GetOptionsDesc(FOptionsDescType Type = OPTIONSDESC_PropertyAssign) const;
   void CopyPropertiesFrom(FWfOptions const &other); // copy data from other object, and emit corresponding signals.
   // kate: syntax c++;
