// this code is GENERATED. Do not change it directly!
public:
   struct FFreeLabelData : public QSharedData {
      QString m_Text;
      FVec3d m_Pos;
      double m_Size;
      uint32_t m_Color;
   };
   typedef QExplicitlySharedDataPointer<FFreeLabelData>
      FSharedDataPtr;
   explicit FFreeLabel(QObject *parent = 0, FSharedDataPtr dref = FSharedDataPtr());
   QSharedPointer<FFreeObject> newLinkedObject(QObject *parent);
   QSharedPointer<FFreeObject> newClonedObject(QObject *parent) const;
protected:
   FSharedDataPtr
      d;
public:
   // update properties in *this from entries in variant map
   Q_SLOT void setp(QVariantMap const &vm);

   Q_PROPERTY(QString text READ GetText WRITE SetText NOTIFY TextChanged)
   Q_SLOT void SetText(QString const &o);
   Q_SIGNAL void TextChanged(QString const &o);
   QString const &GetText() const { return d->m_Text; };

   Q_PROPERTY(FVec3d pos READ GetPos WRITE SetPos NOTIFY PosChanged)
   Q_SLOT void SetPos(FVec3d o);
   Q_SIGNAL void PosChanged(FVec3d const &o);
   FVec3d GetPos() const { return d->m_Pos; };

   Q_PROPERTY(double size READ GetSize WRITE SetSize NOTIFY SizeChanged)
   Q_SLOT void SetSize(double o);
   Q_SIGNAL void SizeChanged(double const &o);
   double GetSize() const { return d->m_Size; };

   Q_PROPERTY(uint32_t color READ GetColor WRITE SetColor NOTIFY ColorChanged)
   Q_SLOT void SetColor(uint32_t o);
   Q_SIGNAL void ColorChanged(uint32_t const &o);
   uint32_t GetColor() const { return d->m_Color; };

   void InitProperties(); // set default values of properties and emit corresponding signals.
   // make script code for setting/collecting options in *this which were changed from their defaults
   QString GetOptionsDesc(FOptionsDescType Type = OPTIONSDESC_PropertyAssign) const;
   void CopyPropertiesFrom(FFreeLabel const &other); // copy data from other object, and emit corresponding signals.
   // kate: syntax c++;
