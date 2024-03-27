// this code is GENERATED. Do not change it directly!
void FElementOptions::InitProperties() {
   m_Color = (int)GetDefaultColor();
   m_BondColor = (int)GetDefaultBondColor();
   m_BondRadiusFactor = 1.3;
   m_CovalentRadius = GetDefaultCovalentRadius();
   m_VdwRadius = GetDefaultVdwRadius();
   m_DrawRadius = GetDefaultDrawRadius();
   m_DrawName = ElementName();
}

QString FElementOptions::GetOptionsDesc(FOptionsDescType Type) const {
   FOptionsDesc desc(Type, "AtomOptions", "doc.elements(ielem)");
   using std::abs;
   if (m_Color != (int)GetDefaultColor())
      desc.setp(-1, "color", m_Color);
   if (m_BondColor != (int)GetDefaultBondColor())
      desc.setp(-1, "bond_color", m_BondColor);
   if (abs(m_BondRadiusFactor - 1.3) > 1e-4)
      desc.setp(-1, "bond_radius_factor", m_BondRadiusFactor);
   if (abs(m_CovalentRadius - GetDefaultCovalentRadius()) > 1e-4)
      desc.setp(-1, "covalent_radius", m_CovalentRadius);
   if (abs(m_VdwRadius - GetDefaultVdwRadius()) > 1e-4)
      desc.setp(-1, "vdw_radius", m_VdwRadius);
   if (abs(m_DrawRadius - GetDefaultDrawRadius()) > 1e-4)
      desc.setp(-1, "draw_radius", m_DrawRadius);
   if (m_DrawName != ElementName())
      desc.setp(-1, "draw_name", m_DrawName);
   return s2q(desc.str());
}

void FElementOptions::setp(QVariantMap const &vm) {
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

void FElementOptions::SetColor(int o) {
   if (o != m_Color) {
      m_Color = o;
      emit ColorChanged(m_Color);
   }
}

void FElementOptions::SetBondColor(int o) {
   if (o != m_BondColor) {
      m_BondColor = o;
      emit BondColorChanged(m_BondColor);
   }
}

void FElementOptions::SetBondRadiusFactor(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_BondRadiusFactor) {
      m_BondRadiusFactor = o;
      emit BondRadiusFactorChanged(m_BondRadiusFactor);
   }
}

void FElementOptions::SetCovalentRadius(double o) {
   if (o < -1.0) o = -1.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_CovalentRadius) {
      m_CovalentRadius = o;
      emit CovalentRadiusChanged(m_CovalentRadius);
   }
}

void FElementOptions::SetVdwRadius(double o) {
   if (o < -1.0) o = -1.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_VdwRadius) {
      m_VdwRadius = o;
      emit VdwRadiusChanged(m_VdwRadius);
   }
}

void FElementOptions::SetDrawRadius(double o) {
   if (o < -1.0) o = -1.0;
   if (o > 4096.0) o = 4096.0;
   if (o != m_DrawRadius) {
      m_DrawRadius = o;
      emit DrawRadiusChanged(m_DrawRadius);
   }
}

void FElementOptions::SetDrawName(QString const &o) {
   if (o != m_DrawName) {
      m_DrawName = o;
      emit DrawNameChanged(m_DrawName);
   }
}

void FElementOptions::CopyPropertiesFrom(FElementOptions const &other) {
   SetColor(other.GetColor());
   SetBondColor(other.GetBondColor());
   SetBondRadiusFactor(other.GetBondRadiusFactor());
   SetCovalentRadius(other.GetCovalentRadius());
   SetVdwRadius(other.GetVdwRadius());
   SetDrawRadius(other.GetDrawRadius());
   SetDrawName(other.GetDrawName());
}
// kate: syntax c++;
