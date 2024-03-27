// this code is GENERATED. Do not change it directly!
FFreeLabel::FFreeLabel(QObject *parent, FSharedDataPtr dref)
  : FFreeObject(parent), d(dref)
{
   if (!this->d) {
      this->d = new FFreeLabelData;
      InitProperties();
   }
}

QSharedPointer<FFreeObject> FFreeLabel::newLinkedObject(QObject *parent) {
   return QSharedPointer<FFreeObject>(new FFreeLabel(parent, this->d));
}

QSharedPointer<FFreeObject> FFreeLabel::newClonedObject(QObject *parent) const {
   return QSharedPointer<FFreeObject>(new FFreeLabel(parent, FSharedDataPtr(new FFreeLabelData(*this->d))));
}

void FFreeLabel::InitProperties() {
   d->m_Text = "0.01";
   d->m_Pos = FVec3d(0,0,0);
   d->m_Size = 0.01;
   d->m_Color = 0xff000000;
}

char const *BoolToCstr(bool o);

QString FFreeLabel::GetOptionsDesc(FOptionsDescType Type) const {
   FOptionsDesc desc(Type, "FreeLabel", "label");
   using std::abs;
   desc.setp(0, "text", d->m_Text);
   desc.setp(1, "pos", d->m_Pos);
   if (abs(d->m_Size - 0.01) > 1e-4)
      desc.setp(-1, "size", d->m_Size);
   if (d->m_Color != 0xff000000)
      desc.setp(-1, "color", d->m_Color);
   return s2q(desc.str());
}

void FFreeLabel::setp(QVariantMap const &vm) {
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

void FFreeLabel::SetText(QString const &o) {
   if (o != d->m_Text) {
      d->m_Text = o;
      emit TextChanged(d->m_Text);
   }
}

void FFreeLabel::SetPos(FVec3d o) {
   if (o != d->m_Pos) {
      d->m_Pos = o;
      emit PosChanged(d->m_Pos);
   }
}

void FFreeLabel::SetSize(double o) {
   if (o < 0) o = 0;
   if (o > 4096) o = 4096;
   if (o != d->m_Size) {
      d->m_Size = o;
      emit SizeChanged(d->m_Size);
   }
}

void FFreeLabel::SetColor(uint32_t o) {
   if (o != d->m_Color) {
      d->m_Color = o;
      emit ColorChanged(d->m_Color);
   }
}

void FFreeLabel::CopyPropertiesFrom(FFreeLabel const &other) {
   SetText(other.GetText());
   SetPos(other.GetPos());
   SetSize(other.GetSize());
   SetColor(other.GetColor());
}
// kate: syntax c++;
