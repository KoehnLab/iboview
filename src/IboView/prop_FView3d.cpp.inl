// this code is GENERATED. Do not change it directly!
void FView3d::InitProperties() {
   m_IsoResolution = 12.0;
   m_IsoThreshold = 80.0;
   m_IsoAutoFlip = false;
   m_AtomScale = 100.0;
   m_BondScale = 100.0;
   m_BondThinning = 0.72;
   m_MultiBondScale = 140.0;
   m_MultiBondPos = 70;
   m_MultiBondType = 1;
   m_BondOrderSource = "Auto";
   m_FullBondThresh = 0.2;
   m_IndicateHiddenAtoms = true;
   m_ShowHydrogens = true;
   m_OrbitalEllipsoidsOnly = false;
   m_OrientOrbitalEllipsoids = true;
   m_FadeType = 1;
   m_FadeWidth = 9.0;
   m_FadeBias = 0.0;
   m_ShaderRegA0 = 0.8;
   m_ShaderRegO0 = 0.8;
   m_ShaderRegA1 = 0.7;
   m_ShaderRegO1 = 0.7;
   m_ShaderRegA2 = 0.4;
   m_ShaderRegO2 = 0.7;
   m_ShaderRegA3 = -0.5;
   m_ShaderRegO3 = -0.5;
   m_ShaderPath = ":/shader";
   m_LabelAtomNumbers = false;
   m_LabelElements = true;
   m_LabelElementsC = true;
   m_LabelElementsH = false;
   m_LabelSize = 80.0;
   m_LabelSizeH = 0.8;
   m_LabelBrightness = -0.2;
   m_DepthPeelingLayers = 4;
   m_SaveAlpha = true;
   m_CropImages = true;
   m_RenderBacksides = false;
   m_SuperSample = true;
   m_FakeAntiAliasing = true;
   connect(this, SIGNAL(IsoResolutionChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(IsoThresholdChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(IsoAutoFlipChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(AtomScaleChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(BondScaleChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(BondThinningChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(MultiBondScaleChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(MultiBondPosChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(MultiBondTypeChanged(int)), this, SLOT(update()));
   connect(this, SIGNAL(BondOrderSourceChanged(QString)), this, SLOT(update()));
   connect(this, SIGNAL(FullBondThreshChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(IndicateHiddenAtomsChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(ShowHydrogensChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(OrbitalEllipsoidsOnlyChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(OrientOrbitalEllipsoidsChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(FadeTypeChanged(int)), this, SLOT(update()));
   connect(this, SIGNAL(FadeWidthChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(FadeBiasChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegA0Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegO0Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegA1Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegO1Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegA2Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegO2Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegA3Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderRegO3Changed(double)), this, SLOT(update()));
   connect(this, SIGNAL(ShaderPathChanged(QString)), this, SLOT(update()));
   connect(this, SIGNAL(LabelAtomNumbersChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(LabelElementsChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(LabelElementsCChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(LabelElementsHChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(LabelSizeChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(LabelSizeHChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(LabelBrightnessChanged(double)), this, SLOT(update()));
   connect(this, SIGNAL(DepthPeelingLayersChanged(int)), this, SLOT(update()));
   connect(this, SIGNAL(SaveAlphaChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(CropImagesChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(RenderBacksidesChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(SuperSampleChanged(bool)), this, SLOT(update()));
   connect(this, SIGNAL(FakeAntiAliasingChanged(bool)), this, SLOT(update()));
}

char const *BoolToCstr(bool o);

QString FView3d::GetOptionsDesc(FOptionsDescType Type) const {
   FOptionsDesc desc(Type, "ViewOptions", "view");
   using std::abs;
   if (abs(m_IsoResolution - 12.0) > 1e-4)
      desc.setp(-1, "iso_resolution", m_IsoResolution);
   if (abs(m_IsoThreshold - 80.0) > 1e-4)
      desc.setp(-1, "iso_threshold", m_IsoThreshold);
   if (m_IsoAutoFlip != false)
      desc.setp(-1, "iso_auto_flip", m_IsoAutoFlip);
   if (abs(m_AtomScale - 100.0) > 1e-4)
      desc.setp(-1, "atom_scale", m_AtomScale);
   if (abs(m_BondScale - 100.0) > 1e-4)
      desc.setp(-1, "bond_scale", m_BondScale);
   if (abs(m_BondThinning - 0.72) > 1e-4)
      desc.setp(-1, "bond_thinning", m_BondThinning);
   if (abs(m_MultiBondScale - 140.0) > 1e-4)
      desc.setp(-1, "multi_bond_scale", m_MultiBondScale);
   if (abs(m_MultiBondPos - 70) > 1e-4)
      desc.setp(-1, "multi_bond_pos", m_MultiBondPos);
   if (m_MultiBondType != 1)
      desc.setp(-1, "multi_bond_type", m_MultiBondType);
   if (m_BondOrderSource != "Auto")
      desc.setp(-1, "bond_order_source", m_BondOrderSource);
   if (abs(m_FullBondThresh - 0.2) > 1e-4)
      desc.setp(-1, "full_bond_thresh", m_FullBondThresh);
   if (m_IndicateHiddenAtoms != true)
      desc.setp(-1, "indicate_hidden_atoms", m_IndicateHiddenAtoms);
   if (m_ShowHydrogens != true)
      desc.setp(-1, "show_hydrogens", m_ShowHydrogens);
   if (m_OrbitalEllipsoidsOnly != false)
      desc.setp(-1, "orbital_ellipsoids_only", m_OrbitalEllipsoidsOnly);
   if (m_OrientOrbitalEllipsoids != true)
      desc.setp(-1, "orient_orbital_ellipsoids", m_OrientOrbitalEllipsoids);
   if (m_FadeType != 1)
      desc.setp(-1, "fade_type", m_FadeType);
   if (abs(m_FadeWidth - 9.0) > 1e-4)
      desc.setp(-1, "fade_width", m_FadeWidth);
   if (abs(m_FadeBias - 0.0) > 1e-4)
      desc.setp(-1, "fade_bias", m_FadeBias);
   if (abs(m_ShaderRegA0 - 0.8) > 1e-4)
      desc.setp(-1, "shader_reg_a0", m_ShaderRegA0);
   if (abs(m_ShaderRegO0 - 0.8) > 1e-4)
      desc.setp(-1, "shader_reg_o0", m_ShaderRegO0);
   if (abs(m_ShaderRegA1 - 0.7) > 1e-4)
      desc.setp(-1, "shader_reg_a1", m_ShaderRegA1);
   if (abs(m_ShaderRegO1 - 0.7) > 1e-4)
      desc.setp(-1, "shader_reg_o1", m_ShaderRegO1);
   if (abs(m_ShaderRegA2 - 0.4) > 1e-4)
      desc.setp(-1, "shader_reg_a2", m_ShaderRegA2);
   if (abs(m_ShaderRegO2 - 0.7) > 1e-4)
      desc.setp(-1, "shader_reg_o2", m_ShaderRegO2);
   if (abs(m_ShaderRegA3 - -0.5) > 1e-4)
      desc.setp(-1, "shader_reg_a3", m_ShaderRegA3);
   if (abs(m_ShaderRegO3 - -0.5) > 1e-4)
      desc.setp(-1, "shader_reg_o3", m_ShaderRegO3);
   if (m_ShaderPath != ":/shader")
      desc.setp(-1, "shader_path", m_ShaderPath);
   if (m_LabelAtomNumbers != false)
      desc.setp(-1, "label_atom_numbers", m_LabelAtomNumbers);
   if (m_LabelElements != true)
      desc.setp(-1, "label_elements", m_LabelElements);
   if (m_LabelElementsC != true)
      desc.setp(-1, "label_elements_c", m_LabelElementsC);
   if (m_LabelElementsH != false)
      desc.setp(-1, "label_elements_h", m_LabelElementsH);
   if (abs(m_LabelSize - 80.0) > 1e-4)
      desc.setp(-1, "label_size", m_LabelSize);
   if (abs(m_LabelSizeH - 0.8) > 1e-4)
      desc.setp(-1, "label_size_h", m_LabelSizeH);
   if (abs(m_LabelBrightness - -0.2) > 1e-4)
      desc.setp(-1, "label_brightness", m_LabelBrightness);
   if (m_DepthPeelingLayers != 4)
      desc.setp(-1, "depth_peeling_layers", m_DepthPeelingLayers);
   if (m_SaveAlpha != true)
      desc.setp(-1, "save_alpha", m_SaveAlpha);
   if (m_CropImages != true)
      desc.setp(-1, "crop_images", m_CropImages);
   if (m_RenderBacksides != false)
      desc.setp(-1, "render_backsides", m_RenderBacksides);
   if (m_SuperSample != true)
      desc.setp(-1, "super_sample", m_SuperSample);
   if (m_FakeAntiAliasing != true)
      desc.setp(-1, "fake_anti_aliasing", m_FakeAntiAliasing);
   return s2q(desc.str());
}

void FView3d::setp(QVariantMap const &vm) {
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

void FView3d::SetIsoResolution(double o) {
   if (o < 1.0) o = 1.0;
   if (o > 100.0) o = 100.0;
   if (o != m_IsoResolution) {
      m_IsoResolution = o;
      emit IsoResolutionChanged(m_IsoResolution);
            IvEmit(" New iso-surfaces will be traced with a resolution of [%1 pt/angstrom]^3.", m_IsoResolution);
   }
}

void FView3d::SetIsoThreshold(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 100.0) o = 100.0;
   if (o != m_IsoThreshold) {
      m_IsoThreshold = o;
      emit IsoThresholdChanged(m_IsoThreshold);
            IvEmit(" New iso-surfaces will be traced with a threshold of %1%.", m_IsoThreshold);
   }
}

void FView3d::SetIsoAutoFlip(bool o) {
   if (o != m_IsoAutoFlip) {
      m_IsoAutoFlip = o;
      emit IsoAutoFlipChanged(m_IsoAutoFlip);
   }
}

void FView3d::SetAtomScale(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_AtomScale) {
      m_AtomScale = o;
      emit AtomScaleChanged(m_AtomScale);
   }
}

void FView3d::SetBondScale(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_BondScale) {
      m_BondScale = o;
      emit BondScaleChanged(m_BondScale);
   }
}

void FView3d::SetBondThinning(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 10.0) o = 10.0;
   if (o != m_BondThinning) {
      m_BondThinning = o;
      emit BondThinningChanged(m_BondThinning);
   }
}

void FView3d::SetMultiBondScale(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_MultiBondScale) {
      m_MultiBondScale = o;
      emit MultiBondScaleChanged(m_MultiBondScale);
   }
}

void FView3d::SetMultiBondPos(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 1000.0) o = 1000.0;
   if (o != m_MultiBondPos) {
      m_MultiBondPos = o;
      emit MultiBondPosChanged(m_MultiBondPos);
   }
}

void FView3d::SetMultiBondType(int o) {
   if (o < 0) o = 0;
   if (o > 1) o = 1;
   if (o != m_MultiBondType) {
      m_MultiBondType = o;
      emit MultiBondTypeChanged(m_MultiBondType);
   }
}

void FView3d::SetBondOrderSource(QString const &o) {
   if (o != m_BondOrderSource) {
      m_BondOrderSource = o;
      emit BondOrderSourceChanged(m_BondOrderSource);
   }
}

void FView3d::SetFullBondThresh(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 0.5) o = 0.5;
   if (o != m_FullBondThresh) {
      m_FullBondThresh = o;
      emit FullBondThreshChanged(m_FullBondThresh);
   }
}

void FView3d::SetIndicateHiddenAtoms(bool o) {
   if (o != m_IndicateHiddenAtoms) {
      m_IndicateHiddenAtoms = o;
      emit IndicateHiddenAtomsChanged(m_IndicateHiddenAtoms);
   }
}

void FView3d::SetShowHydrogens(bool o) {
   if (o != m_ShowHydrogens) {
      m_ShowHydrogens = o;
      emit ShowHydrogensChanged(m_ShowHydrogens);
   }
}

void FView3d::SetOrbitalEllipsoidsOnly(bool o) {
   if (o != m_OrbitalEllipsoidsOnly) {
      m_OrbitalEllipsoidsOnly = o;
      emit OrbitalEllipsoidsOnlyChanged(m_OrbitalEllipsoidsOnly);
   }
}

void FView3d::SetOrientOrbitalEllipsoids(bool o) {
   if (o != m_OrientOrbitalEllipsoids) {
      m_OrientOrbitalEllipsoids = o;
      emit OrientOrbitalEllipsoidsChanged(m_OrientOrbitalEllipsoids);
   }
}

void FView3d::SetFadeType(int o) {
   if (o < 0) o = 0;
   if (o > 1) o = 1;
   if (o != m_FadeType) {
      m_FadeType = o;
      emit FadeTypeChanged(m_FadeType);
   }
}

void FView3d::SetFadeWidth(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_FadeWidth) {
      m_FadeWidth = o;
      emit FadeWidthChanged(m_FadeWidth);
   }
}

void FView3d::SetFadeBias(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_FadeBias) {
      m_FadeBias = o;
      emit FadeBiasChanged(m_FadeBias);
   }
}

void FView3d::SetShaderRegA0(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegA0) {
      m_ShaderRegA0 = o;
      emit ShaderRegA0Changed(m_ShaderRegA0);
   }
}

void FView3d::SetShaderRegO0(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegO0) {
      m_ShaderRegO0 = o;
      emit ShaderRegO0Changed(m_ShaderRegO0);
   }
}

void FView3d::SetShaderRegA1(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegA1) {
      m_ShaderRegA1 = o;
      emit ShaderRegA1Changed(m_ShaderRegA1);
   }
}

void FView3d::SetShaderRegO1(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegO1) {
      m_ShaderRegO1 = o;
      emit ShaderRegO1Changed(m_ShaderRegO1);
   }
}

void FView3d::SetShaderRegA2(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegA2) {
      m_ShaderRegA2 = o;
      emit ShaderRegA2Changed(m_ShaderRegA2);
   }
}

void FView3d::SetShaderRegO2(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegO2) {
      m_ShaderRegO2 = o;
      emit ShaderRegO2Changed(m_ShaderRegO2);
   }
}

void FView3d::SetShaderRegA3(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegA3) {
      m_ShaderRegA3 = o;
      emit ShaderRegA3Changed(m_ShaderRegA3);
   }
}

void FView3d::SetShaderRegO3(double o) {
   if (o < -9999.0) o = -9999.0;
   if (o > 9999.0) o = 9999.0;
   if (o != m_ShaderRegO3) {
      m_ShaderRegO3 = o;
      emit ShaderRegO3Changed(m_ShaderRegO3);
   }
}

void FView3d::SetShaderPath(QString const &o) {
   if (o != m_ShaderPath) {
      m_ShaderPath = o;
      emit ShaderPathChanged(m_ShaderPath);
   }
}

void FView3d::SetLabelAtomNumbers(bool o) {
   if (o != m_LabelAtomNumbers) {
      m_LabelAtomNumbers = o;
      emit LabelAtomNumbersChanged(m_LabelAtomNumbers);
   }
}

void FView3d::SetLabelElements(bool o) {
   if (o != m_LabelElements) {
      m_LabelElements = o;
      emit LabelElementsChanged(m_LabelElements);
   }
}

void FView3d::SetLabelElementsC(bool o) {
   if (o != m_LabelElementsC) {
      m_LabelElementsC = o;
      emit LabelElementsCChanged(m_LabelElementsC);
   }
}

void FView3d::SetLabelElementsH(bool o) {
   if (o != m_LabelElementsH) {
      m_LabelElementsH = o;
      emit LabelElementsHChanged(m_LabelElementsH);
   }
}

void FView3d::SetLabelSize(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 1000.0) o = 1000.0;
   if (o != m_LabelSize) {
      m_LabelSize = o;
      emit LabelSizeChanged(m_LabelSize);
   }
}

void FView3d::SetLabelSizeH(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 10.0) o = 10.0;
   if (o != m_LabelSizeH) {
      m_LabelSizeH = o;
      emit LabelSizeHChanged(m_LabelSizeH);
   }
}

void FView3d::SetLabelBrightness(double o) {
   if (o < -1.0) o = -1.0;
   if (o > 1.0) o = 1.0;
   if (o != m_LabelBrightness) {
      m_LabelBrightness = o;
      emit LabelBrightnessChanged(m_LabelBrightness);
   }
}

void FView3d::SetDepthPeelingLayers(int o) {
   if (o < 0) o = 0;
   if (o > 16) o = 16;
   if (o != m_DepthPeelingLayers) {
      m_DepthPeelingLayers = o;
      emit DepthPeelingLayersChanged(m_DepthPeelingLayers);
   }
}

void FView3d::SetSaveAlpha(bool o) {
   if (o != m_SaveAlpha) {
      m_SaveAlpha = o;
      emit SaveAlphaChanged(m_SaveAlpha);
   }
}

void FView3d::SetCropImages(bool o) {
   if (o != m_CropImages) {
      m_CropImages = o;
      emit CropImagesChanged(m_CropImages);
   }
}

void FView3d::SetRenderBacksides(bool o) {
   if (o != m_RenderBacksides) {
      m_RenderBacksides = o;
      emit RenderBacksidesChanged(m_RenderBacksides);
   }
}

void FView3d::SetSuperSample(bool o) {
   if (o != m_SuperSample) {
      m_SuperSample = o;
      emit SuperSampleChanged(m_SuperSample);
   }
}

void FView3d::SetFakeAntiAliasing(bool o) {
   if (o != m_FakeAntiAliasing) {
      m_FakeAntiAliasing = o;
      emit FakeAntiAliasingChanged(m_FakeAntiAliasing);
   }
}

void FView3d::CopyPropertiesFrom(FView3d const &other) {
   SetIsoResolution(other.GetIsoResolution());
   SetIsoThreshold(other.GetIsoThreshold());
   SetIsoAutoFlip(other.GetIsoAutoFlip());
   SetAtomScale(other.GetAtomScale());
   SetBondScale(other.GetBondScale());
   SetBondThinning(other.GetBondThinning());
   SetMultiBondScale(other.GetMultiBondScale());
   SetMultiBondPos(other.GetMultiBondPos());
   SetMultiBondType(other.GetMultiBondType());
   SetBondOrderSource(other.GetBondOrderSource());
   SetFullBondThresh(other.GetFullBondThresh());
   SetIndicateHiddenAtoms(other.GetIndicateHiddenAtoms());
   SetShowHydrogens(other.GetShowHydrogens());
   SetOrbitalEllipsoidsOnly(other.GetOrbitalEllipsoidsOnly());
   SetOrientOrbitalEllipsoids(other.GetOrientOrbitalEllipsoids());
   SetFadeType(other.GetFadeType());
   SetFadeWidth(other.GetFadeWidth());
   SetFadeBias(other.GetFadeBias());
   SetShaderRegA0(other.GetShaderRegA0());
   SetShaderRegO0(other.GetShaderRegO0());
   SetShaderRegA1(other.GetShaderRegA1());
   SetShaderRegO1(other.GetShaderRegO1());
   SetShaderRegA2(other.GetShaderRegA2());
   SetShaderRegO2(other.GetShaderRegO2());
   SetShaderRegA3(other.GetShaderRegA3());
   SetShaderRegO3(other.GetShaderRegO3());
   SetShaderPath(other.GetShaderPath());
   SetLabelAtomNumbers(other.GetLabelAtomNumbers());
   SetLabelElements(other.GetLabelElements());
   SetLabelElementsC(other.GetLabelElementsC());
   SetLabelElementsH(other.GetLabelElementsH());
   SetLabelSize(other.GetLabelSize());
   SetLabelSizeH(other.GetLabelSizeH());
   SetLabelBrightness(other.GetLabelBrightness());
   SetDepthPeelingLayers(other.GetDepthPeelingLayers());
   SetSaveAlpha(other.GetSaveAlpha());
   SetCropImages(other.GetCropImages());
   SetRenderBacksides(other.GetRenderBacksides());
   SetSuperSample(other.GetSuperSample());
   SetFakeAntiAliasing(other.GetFakeAntiAliasing());
}
// kate: syntax c++;
