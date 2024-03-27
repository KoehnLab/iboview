// this code is GENERATED. Do not change it directly!
   // update properties in *this from entries in variant map
   Q_SLOT void setp(QVariantMap const &vm);

   Q_PROPERTY(int color READ GetColor WRITE SetColor NOTIFY ColorChanged)
   int m_Color;
   Q_SLOT void SetColor(int o);
   Q_SIGNAL void ColorChanged(int const &o);
   int GetColor() const { return m_Color; };

   Q_PROPERTY(int bond_color READ GetBondColor WRITE SetBondColor NOTIFY BondColorChanged)
   int m_BondColor;
   Q_SLOT void SetBondColor(int o);
   Q_SIGNAL void BondColorChanged(int const &o);
   int GetBondColor() const { return m_BondColor; };

   Q_PROPERTY(double bond_radius_factor READ GetBondRadiusFactor WRITE SetBondRadiusFactor NOTIFY BondRadiusFactorChanged)
   double m_BondRadiusFactor;
   Q_SLOT void SetBondRadiusFactor(double o);
   Q_SIGNAL void BondRadiusFactorChanged(double const &o);
   double GetBondRadiusFactor() const { return m_BondRadiusFactor; };

   Q_PROPERTY(double covalent_radius READ GetCovalentRadius WRITE SetCovalentRadius NOTIFY CovalentRadiusChanged)
   double m_CovalentRadius;
   Q_SLOT void SetCovalentRadius(double o);
   Q_SIGNAL void CovalentRadiusChanged(double const &o);
   double GetCovalentRadius() const { return m_CovalentRadius; };

   Q_PROPERTY(double vdw_radius READ GetVdwRadius WRITE SetVdwRadius NOTIFY VdwRadiusChanged)
   double m_VdwRadius;
   Q_SLOT void SetVdwRadius(double o);
   Q_SIGNAL void VdwRadiusChanged(double const &o);
   double GetVdwRadius() const { return m_VdwRadius; };

   Q_PROPERTY(double draw_radius READ GetDrawRadius WRITE SetDrawRadius NOTIFY DrawRadiusChanged)
   double m_DrawRadius;
   Q_SLOT void SetDrawRadius(double o);
   Q_SIGNAL void DrawRadiusChanged(double const &o);
   double GetDrawRadius() const { return m_DrawRadius; };

   Q_PROPERTY(QString draw_name READ GetDrawName WRITE SetDrawName NOTIFY DrawNameChanged)
   QString m_DrawName;
   Q_SLOT void SetDrawName(QString const &o);
   Q_SIGNAL void DrawNameChanged(QString const &o);
   QString const &GetDrawName() const { return m_DrawName; };

   void InitProperties(); // set default values of properties and emit corresponding signals.
   // make script code for setting/collecting options in *this which were changed from their defaults
   QString GetOptionsDesc(FOptionsDescType Type = OPTIONSDESC_PropertyAssign) const;
   void CopyPropertiesFrom(FElementOptions const &other); // copy data from other object, and emit corresponding signals.
   // kate: syntax c++;
