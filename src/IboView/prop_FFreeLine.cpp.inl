// this code is GENERATED. Do not change it directly!
FFreeLine::FFreeLine(QObject *parent, FSharedDataPtr dref)
  : FFreeObject(parent), d(dref)
{
   if (!this->d) {
      this->d = new FFreeLineData;
      InitProperties();
   }
}

QSharedPointer<FFreeObject> FFreeLine::newLinkedObject(QObject *parent) {
   return QSharedPointer<FFreeObject>(new FFreeLine(parent, this->d));
}

QSharedPointer<FFreeObject> FFreeLine::newClonedObject(QObject *parent) const {
   return QSharedPointer<FFreeObject>(new FFreeLine(parent, FSharedDataPtr(new FFreeLineData(*this->d))));
}

void FFreeLine::InitProperties() {
   d->m_From = FVec3d(0,0,0);
   d->m_To = FVec3d(0,0,0);
   d->m_Color = 0xff000000;
   d->m_Width = 0.01;
   d->m_Weight = 1.0;
}

char const *BoolToCstr(bool o);

QString FFreeLine::GetOptionsDesc(FOptionsDescType Type) const {
   FOptionsDesc desc(Type, "FreeLine", "line");
   using std::abs;
   desc.setp(0, "from", d->m_From);
   desc.setp(1, "to", d->m_To);
   if (d->m_Color != 0xff000000)
      desc.setp(-1, "color", d->m_Color);
   if (abs(d->m_Width - 0.01) > 1e-4)
      desc.setp(-1, "width", d->m_Width);
   if (abs(d->m_Weight - 1.0) > 1e-4)
      desc.setp(-1, "weight", d->m_Weight);
   return s2q(desc.str());
}

void FFreeLine::setp(QVariantMap const &vm) {
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

void FFreeLine::SetFrom(FVec3d o) {
   if (o != d->m_From) {
      d->m_From = o;
      emit FromChanged(d->m_From);
   }
}

void FFreeLine::SetTo(FVec3d o) {
   if (o != d->m_To) {
      d->m_To = o;
      emit ToChanged(d->m_To);
   }
}

void FFreeLine::SetColor(uint32_t o) {
   if (o != d->m_Color) {
      d->m_Color = o;
      emit ColorChanged(d->m_Color);
   }
}

void FFreeLine::SetWidth(double o) {
   if (o < 0) o = 0;
   if (o > 4096) o = 4096;
   if (o != d->m_Width) {
      d->m_Width = o;
      emit WidthChanged(d->m_Width);
   }
}

void FFreeLine::SetWeight(double o) {
   if (o < 0.0) o = 0.0;
   if (o > 1.0) o = 1.0;
   if (o != d->m_Weight) {
      d->m_Weight = o;
      emit WeightChanged(d->m_Weight);
   }
}

void FFreeLine::CopyPropertiesFrom(FFreeLine const &other) {
   SetFrom(other.GetFrom());
   SetTo(other.GetTo());
   SetColor(other.GetColor());
   SetWidth(other.GetWidth());
   SetWeight(other.GetWeight());
}
// kate: syntax c++;
