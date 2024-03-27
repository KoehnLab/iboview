// this code is GENERATED. Do not change it directly!
   // update properties in *this from entries in variant map
   Q_SLOT void setp(QVariantMap const &vm);

   Q_PROPERTY(double iso_resolution READ GetIsoResolution WRITE SetIsoResolution NOTIFY IsoResolutionChanged)
   double m_IsoResolution;
   Q_SLOT void SetIsoResolution(double o);
   Q_SIGNAL void IsoResolutionChanged(double const &o);
   double GetIsoResolution() const { return m_IsoResolution; };

   Q_PROPERTY(double iso_threshold READ GetIsoThreshold WRITE SetIsoThreshold NOTIFY IsoThresholdChanged)
   double m_IsoThreshold;
   Q_SLOT void SetIsoThreshold(double o);
   Q_SIGNAL void IsoThresholdChanged(double const &o);
   double GetIsoThreshold() const { return m_IsoThreshold; };

   Q_PROPERTY(bool iso_auto_flip READ GetIsoAutoFlip WRITE SetIsoAutoFlip NOTIFY IsoAutoFlipChanged)
   bool m_IsoAutoFlip;
   Q_SLOT void SetIsoAutoFlip(bool o);
   Q_SIGNAL void IsoAutoFlipChanged(bool const &o);
   bool GetIsoAutoFlip() const { return m_IsoAutoFlip; };

   Q_PROPERTY(double atom_scale READ GetAtomScale WRITE SetAtomScale NOTIFY AtomScaleChanged)
   double m_AtomScale;
   Q_SLOT void SetAtomScale(double o);
   Q_SIGNAL void AtomScaleChanged(double const &o);
   double GetAtomScale() const { return m_AtomScale; };

   Q_PROPERTY(double bond_scale READ GetBondScale WRITE SetBondScale NOTIFY BondScaleChanged)
   double m_BondScale;
   Q_SLOT void SetBondScale(double o);
   Q_SIGNAL void BondScaleChanged(double const &o);
   double GetBondScale() const { return m_BondScale; };

   Q_PROPERTY(double bond_thinning READ GetBondThinning WRITE SetBondThinning NOTIFY BondThinningChanged)
   double m_BondThinning;
   Q_SLOT void SetBondThinning(double o);
   Q_SIGNAL void BondThinningChanged(double const &o);
   double GetBondThinning() const { return m_BondThinning; };

   Q_PROPERTY(double multi_bond_scale READ GetMultiBondScale WRITE SetMultiBondScale NOTIFY MultiBondScaleChanged)
   double m_MultiBondScale;
   Q_SLOT void SetMultiBondScale(double o);
   Q_SIGNAL void MultiBondScaleChanged(double const &o);
   double GetMultiBondScale() const { return m_MultiBondScale; };

   Q_PROPERTY(double multi_bond_pos READ GetMultiBondPos WRITE SetMultiBondPos NOTIFY MultiBondPosChanged)
   double m_MultiBondPos;
   Q_SLOT void SetMultiBondPos(double o);
   Q_SIGNAL void MultiBondPosChanged(double const &o);
   double GetMultiBondPos() const { return m_MultiBondPos; };

   Q_PROPERTY(int multi_bond_type READ GetMultiBondType WRITE SetMultiBondType NOTIFY MultiBondTypeChanged)
   int m_MultiBondType;
   Q_SLOT void SetMultiBondType(int o);
   Q_SIGNAL void MultiBondTypeChanged(int const &o);
   int GetMultiBondType() const { return m_MultiBondType; };

   Q_PROPERTY(QString bond_order_source READ GetBondOrderSource WRITE SetBondOrderSource NOTIFY BondOrderSourceChanged)
   QString m_BondOrderSource;
   Q_SLOT void SetBondOrderSource(QString const &o);
   Q_SIGNAL void BondOrderSourceChanged(QString const &o);
   QString const &GetBondOrderSource() const { return m_BondOrderSource; };

   Q_PROPERTY(double full_bond_thresh READ GetFullBondThresh WRITE SetFullBondThresh NOTIFY FullBondThreshChanged)
   double m_FullBondThresh;
   Q_SLOT void SetFullBondThresh(double o);
   Q_SIGNAL void FullBondThreshChanged(double const &o);
   double GetFullBondThresh() const { return m_FullBondThresh; };

   Q_PROPERTY(bool indicate_hidden_atoms READ GetIndicateHiddenAtoms WRITE SetIndicateHiddenAtoms NOTIFY IndicateHiddenAtomsChanged)
   bool m_IndicateHiddenAtoms;
   Q_SLOT void SetIndicateHiddenAtoms(bool o);
   Q_SIGNAL void IndicateHiddenAtomsChanged(bool const &o);
   bool GetIndicateHiddenAtoms() const { return m_IndicateHiddenAtoms; };

   Q_PROPERTY(bool show_hydrogens READ GetShowHydrogens WRITE SetShowHydrogens NOTIFY ShowHydrogensChanged)
   bool m_ShowHydrogens;
   Q_SLOT void SetShowHydrogens(bool o);
   Q_SIGNAL void ShowHydrogensChanged(bool const &o);
   bool GetShowHydrogens() const { return m_ShowHydrogens; };

   Q_PROPERTY(bool orbital_ellipsoids_only READ GetOrbitalEllipsoidsOnly WRITE SetOrbitalEllipsoidsOnly NOTIFY OrbitalEllipsoidsOnlyChanged)
   bool m_OrbitalEllipsoidsOnly;
   Q_SLOT void SetOrbitalEllipsoidsOnly(bool o);
   Q_SIGNAL void OrbitalEllipsoidsOnlyChanged(bool const &o);
   bool GetOrbitalEllipsoidsOnly() const { return m_OrbitalEllipsoidsOnly; };

   Q_PROPERTY(bool orient_orbital_ellipsoids READ GetOrientOrbitalEllipsoids WRITE SetOrientOrbitalEllipsoids NOTIFY OrientOrbitalEllipsoidsChanged)
   bool m_OrientOrbitalEllipsoids;
   Q_SLOT void SetOrientOrbitalEllipsoids(bool o);
   Q_SIGNAL void OrientOrbitalEllipsoidsChanged(bool const &o);
   bool GetOrientOrbitalEllipsoids() const { return m_OrientOrbitalEllipsoids; };

   Q_PROPERTY(int fade_type READ GetFadeType WRITE SetFadeType NOTIFY FadeTypeChanged)
   int m_FadeType;
   Q_SLOT void SetFadeType(int o);
   Q_SIGNAL void FadeTypeChanged(int const &o);
   int GetFadeType() const { return m_FadeType; };

   Q_PROPERTY(double fade_width READ GetFadeWidth WRITE SetFadeWidth NOTIFY FadeWidthChanged)
   double m_FadeWidth;
   Q_SLOT void SetFadeWidth(double o);
   Q_SIGNAL void FadeWidthChanged(double const &o);
   double GetFadeWidth() const { return m_FadeWidth; };

   Q_PROPERTY(double fade_bias READ GetFadeBias WRITE SetFadeBias NOTIFY FadeBiasChanged)
   double m_FadeBias;
   Q_SLOT void SetFadeBias(double o);
   Q_SIGNAL void FadeBiasChanged(double const &o);
   double GetFadeBias() const { return m_FadeBias; };

   Q_PROPERTY(double shader_reg_a0 READ GetShaderRegA0 WRITE SetShaderRegA0 NOTIFY ShaderRegA0Changed)
   double m_ShaderRegA0;
   Q_SLOT void SetShaderRegA0(double o);
   Q_SIGNAL void ShaderRegA0Changed(double const &o);
   double GetShaderRegA0() const { return m_ShaderRegA0; };

   Q_PROPERTY(double shader_reg_o0 READ GetShaderRegO0 WRITE SetShaderRegO0 NOTIFY ShaderRegO0Changed)
   double m_ShaderRegO0;
   Q_SLOT void SetShaderRegO0(double o);
   Q_SIGNAL void ShaderRegO0Changed(double const &o);
   double GetShaderRegO0() const { return m_ShaderRegO0; };

   Q_PROPERTY(double shader_reg_a1 READ GetShaderRegA1 WRITE SetShaderRegA1 NOTIFY ShaderRegA1Changed)
   double m_ShaderRegA1;
   Q_SLOT void SetShaderRegA1(double o);
   Q_SIGNAL void ShaderRegA1Changed(double const &o);
   double GetShaderRegA1() const { return m_ShaderRegA1; };

   Q_PROPERTY(double shader_reg_o1 READ GetShaderRegO1 WRITE SetShaderRegO1 NOTIFY ShaderRegO1Changed)
   double m_ShaderRegO1;
   Q_SLOT void SetShaderRegO1(double o);
   Q_SIGNAL void ShaderRegO1Changed(double const &o);
   double GetShaderRegO1() const { return m_ShaderRegO1; };

   Q_PROPERTY(double shader_reg_a2 READ GetShaderRegA2 WRITE SetShaderRegA2 NOTIFY ShaderRegA2Changed)
   double m_ShaderRegA2;
   Q_SLOT void SetShaderRegA2(double o);
   Q_SIGNAL void ShaderRegA2Changed(double const &o);
   double GetShaderRegA2() const { return m_ShaderRegA2; };

   Q_PROPERTY(double shader_reg_o2 READ GetShaderRegO2 WRITE SetShaderRegO2 NOTIFY ShaderRegO2Changed)
   double m_ShaderRegO2;
   Q_SLOT void SetShaderRegO2(double o);
   Q_SIGNAL void ShaderRegO2Changed(double const &o);
   double GetShaderRegO2() const { return m_ShaderRegO2; };

   Q_PROPERTY(double shader_reg_a3 READ GetShaderRegA3 WRITE SetShaderRegA3 NOTIFY ShaderRegA3Changed)
   double m_ShaderRegA3;
   Q_SLOT void SetShaderRegA3(double o);
   Q_SIGNAL void ShaderRegA3Changed(double const &o);
   double GetShaderRegA3() const { return m_ShaderRegA3; };

   Q_PROPERTY(double shader_reg_o3 READ GetShaderRegO3 WRITE SetShaderRegO3 NOTIFY ShaderRegO3Changed)
   double m_ShaderRegO3;
   Q_SLOT void SetShaderRegO3(double o);
   Q_SIGNAL void ShaderRegO3Changed(double const &o);
   double GetShaderRegO3() const { return m_ShaderRegO3; };

   Q_PROPERTY(QString shader_path READ GetShaderPath WRITE SetShaderPath NOTIFY ShaderPathChanged)
   QString m_ShaderPath;
   Q_SLOT void SetShaderPath(QString const &o);
   Q_SIGNAL void ShaderPathChanged(QString const &o);
   QString const &GetShaderPath() const { return m_ShaderPath; };

   Q_PROPERTY(bool label_atom_numbers READ GetLabelAtomNumbers WRITE SetLabelAtomNumbers NOTIFY LabelAtomNumbersChanged)
   bool m_LabelAtomNumbers;
   Q_SLOT void SetLabelAtomNumbers(bool o);
   Q_SIGNAL void LabelAtomNumbersChanged(bool const &o);
   bool GetLabelAtomNumbers() const { return m_LabelAtomNumbers; };

   Q_PROPERTY(bool label_elements READ GetLabelElements WRITE SetLabelElements NOTIFY LabelElementsChanged)
   bool m_LabelElements;
   Q_SLOT void SetLabelElements(bool o);
   Q_SIGNAL void LabelElementsChanged(bool const &o);
   bool GetLabelElements() const { return m_LabelElements; };

   Q_PROPERTY(bool label_elements_c READ GetLabelElementsC WRITE SetLabelElementsC NOTIFY LabelElementsCChanged)
   bool m_LabelElementsC;
   Q_SLOT void SetLabelElementsC(bool o);
   Q_SIGNAL void LabelElementsCChanged(bool const &o);
   bool GetLabelElementsC() const { return m_LabelElementsC; };

   Q_PROPERTY(bool label_elements_h READ GetLabelElementsH WRITE SetLabelElementsH NOTIFY LabelElementsHChanged)
   bool m_LabelElementsH;
   Q_SLOT void SetLabelElementsH(bool o);
   Q_SIGNAL void LabelElementsHChanged(bool const &o);
   bool GetLabelElementsH() const { return m_LabelElementsH; };

   Q_PROPERTY(double label_size READ GetLabelSize WRITE SetLabelSize NOTIFY LabelSizeChanged)
   double m_LabelSize;
   Q_SLOT void SetLabelSize(double o);
   Q_SIGNAL void LabelSizeChanged(double const &o);
   double GetLabelSize() const { return m_LabelSize; };

   Q_PROPERTY(double label_size_h READ GetLabelSizeH WRITE SetLabelSizeH NOTIFY LabelSizeHChanged)
   double m_LabelSizeH;
   Q_SLOT void SetLabelSizeH(double o);
   Q_SIGNAL void LabelSizeHChanged(double const &o);
   double GetLabelSizeH() const { return m_LabelSizeH; };

   Q_PROPERTY(double label_brightness READ GetLabelBrightness WRITE SetLabelBrightness NOTIFY LabelBrightnessChanged)
   double m_LabelBrightness;
   Q_SLOT void SetLabelBrightness(double o);
   Q_SIGNAL void LabelBrightnessChanged(double const &o);
   double GetLabelBrightness() const { return m_LabelBrightness; };

   Q_PROPERTY(int depth_peeling_layers READ GetDepthPeelingLayers WRITE SetDepthPeelingLayers NOTIFY DepthPeelingLayersChanged)
   int m_DepthPeelingLayers;
   Q_SLOT void SetDepthPeelingLayers(int o);
   Q_SIGNAL void DepthPeelingLayersChanged(int const &o);
   int GetDepthPeelingLayers() const { return m_DepthPeelingLayers; };

   Q_PROPERTY(bool save_alpha READ GetSaveAlpha WRITE SetSaveAlpha NOTIFY SaveAlphaChanged)
   bool m_SaveAlpha;
   Q_SLOT void SetSaveAlpha(bool o);
   Q_SIGNAL void SaveAlphaChanged(bool const &o);
   bool GetSaveAlpha() const { return m_SaveAlpha; };

   Q_PROPERTY(bool crop_images READ GetCropImages WRITE SetCropImages NOTIFY CropImagesChanged)
   bool m_CropImages;
   Q_SLOT void SetCropImages(bool o);
   Q_SIGNAL void CropImagesChanged(bool const &o);
   bool GetCropImages() const { return m_CropImages; };

   Q_PROPERTY(bool render_backsides READ GetRenderBacksides WRITE SetRenderBacksides NOTIFY RenderBacksidesChanged)
   bool m_RenderBacksides;
   Q_SLOT void SetRenderBacksides(bool o);
   Q_SIGNAL void RenderBacksidesChanged(bool const &o);
   bool GetRenderBacksides() const { return m_RenderBacksides; };

   Q_PROPERTY(bool super_sample READ GetSuperSample WRITE SetSuperSample NOTIFY SuperSampleChanged)
   bool m_SuperSample;
   Q_SLOT void SetSuperSample(bool o);
   Q_SIGNAL void SuperSampleChanged(bool const &o);
   bool GetSuperSample() const { return m_SuperSample; };

   Q_PROPERTY(bool fake_anti_aliasing READ GetFakeAntiAliasing WRITE SetFakeAntiAliasing NOTIFY FakeAntiAliasingChanged)
   bool m_FakeAntiAliasing;
   Q_SLOT void SetFakeAntiAliasing(bool o);
   Q_SIGNAL void FakeAntiAliasingChanged(bool const &o);
   bool GetFakeAntiAliasing() const { return m_FakeAntiAliasing; };

   void InitProperties(); // set default values of properties and emit corresponding signals.
   // make script code for setting/collecting options in *this which were changed from their defaults
   QString GetOptionsDesc(FOptionsDescType Type = OPTIONSDESC_PropertyAssign) const;
   void CopyPropertiesFrom(FView3d const &other); // copy data from other object, and emit corresponding signals.
   // kate: syntax c++;
