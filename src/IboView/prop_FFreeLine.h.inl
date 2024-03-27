// this code is GENERATED. Do not change it directly!
public:
   struct FFreeLineData : public QSharedData {
      FVec3d m_From;
      FVec3d m_To;
      uint32_t m_Color;
      double m_Width;
      double m_Weight;
   };
   typedef QExplicitlySharedDataPointer<FFreeLineData>
      FSharedDataPtr;
   explicit FFreeLine(QObject *parent = 0, FSharedDataPtr dref = FSharedDataPtr());
   QSharedPointer<FFreeObject> newLinkedObject(QObject *parent);
   QSharedPointer<FFreeObject> newClonedObject(QObject *parent) const;
protected:
   FSharedDataPtr
      d;
public:
   // update properties in *this from entries in variant map
   Q_SLOT void setp(QVariantMap const &vm);

   Q_PROPERTY(FVec3d from READ GetFrom WRITE SetFrom NOTIFY FromChanged)
   Q_SLOT void SetFrom(FVec3d o);
   Q_SIGNAL void FromChanged(FVec3d const &o);
   FVec3d GetFrom() const { return d->m_From; };

   Q_PROPERTY(FVec3d to READ GetTo WRITE SetTo NOTIFY ToChanged)
   Q_SLOT void SetTo(FVec3d o);
   Q_SIGNAL void ToChanged(FVec3d const &o);
   FVec3d GetTo() const { return d->m_To; };

   Q_PROPERTY(uint32_t color READ GetColor WRITE SetColor NOTIFY ColorChanged)
   Q_SLOT void SetColor(uint32_t o);
   Q_SIGNAL void ColorChanged(uint32_t const &o);
   uint32_t GetColor() const { return d->m_Color; };

   Q_PROPERTY(double width READ GetWidth WRITE SetWidth NOTIFY WidthChanged)
   Q_SLOT void SetWidth(double o);
   Q_SIGNAL void WidthChanged(double const &o);
   double GetWidth() const { return d->m_Width; };

   Q_PROPERTY(double weight READ GetWeight WRITE SetWeight NOTIFY WeightChanged)
   Q_SLOT void SetWeight(double o);
   Q_SIGNAL void WeightChanged(double const &o);
   double GetWeight() const { return d->m_Weight; };

   void InitProperties(); // set default values of properties and emit corresponding signals.
   // make script code for setting/collecting options in *this which were changed from their defaults
   QString GetOptionsDesc(FOptionsDescType Type = OPTIONSDESC_PropertyAssign) const;
   void CopyPropertiesFrom(FFreeLine const &other); // copy data from other object, and emit corresponding signals.
   // kate: syntax c++;
