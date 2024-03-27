/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#ifndef ORBVIEW_MAIN_H
#define ORBVIEW_MAIN_H

#include "Iv.h"
#include <QMainWindow>
#include <QAction>
#include <QSortFilterProxyModel>
#include <QScriptEngine>
#include <QResizeEvent>
#include <QShowEvent>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QCloseEvent>
#include <QDialog>
#include <QEvent>
#include "IvView3D.h"
#include "IvDocument.h"
// #include "IvScript.h"

namespace Ui{
    class MainWindow;
    class AboutDialog;
}

class FDocument;
class IApplication;
class IView3d;

class FShowOnlyCurrentColumnFilter : public QSortFilterProxyModel
{
   Q_OBJECT;
public:
   typedef QSortFilterProxyModel
      FBase;
   bool filterAcceptsColumn(int sourceColumn, const QModelIndex &sourceParent) const; // override
   explicit FShowOnlyCurrentColumnFilter(int ShownColumn_, QObject *parent = 0);

   int shownColumn() const { return m_ShownColumn; }
public slots:
   void setShownColumn(int NewColumn);
private:
//    FDocument *m_pDocument;
   int
      m_ShownColumn;
};


// class FMainWindow : public QMainWindow, public IApplication
class FMainWindow : public QMainWindow
{
    Q_OBJECT
public:
   typedef QMainWindow
      FBase;
   FMainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
   ~FMainWindow();
public:
   Ui::MainWindow *ui;
   IView3d *view3d;
   FDocument *document;
   FShowOnlyCurrentColumnFilter *dataViewFilter;

public:
   FGeometry *pGetActiveGeometry();
   FOrbital *pGetOrbital(int iMo);

   QString MakeStateScript();
   void WriteStateScript(QString FileName="");

   // fixme: remove this if we get a saturation control for main color.
//    float hc, vc, sc, ac;
   int iUpdateLocked;

public slots:
   void onMoveFrameButtonClicked();
   void onDataChanged();
   void onActiveDatasetChanged();
   void onFrameIdChanged(int iNewFrame);
   void onDataRowClicked(QModelIndex const &);
   void onDataRowDoubleClicked(QModelIndex const &);
   void onActiveIboCurveChanged();

// various more-or-less directly UI-exposed actions.
   void onOpenFile();
   void onSaveState();
   void onSaveStateAs();
   void onCopyState();
   void onExecStateFromClipboard();
   void onSavePicture();
   void onSavePictureAs();
   void onCopyPicture();
   void onExportCurves();
   void onOpenPreferences();
   void onFindReactingOrbitals();

   void onTraceIsoSurfacesClicked();
   void onAboutClicked();
   void onToggleTrackedOrbitalClicked();
   void onComputeWaveFunctionTriggered();
   void ExecPresetScript();
   void toggleViewPropertyViaAction(bool newState);

   void onRebuildIsoSurfacesClicked();
   void onFlipOrbitalClicked();
   void toggleShadingControls();
   void setShaderPath();

   void onOrbitalColorChanged(int iNewValue);
   void onDiscreeteTrafoTriggered();
   void processInputs();

   void showFrameLog();
   void onTitleChanged();
   void UpdateAtomAndBondScalesText();
   void ShowTablesAndCurves();
   void ShowEditFramesForm();
   void ShowEditVolumeSurfacesForm();
   void ShowComputeEosForm();

   void dummySignalTest(double o);

   void StoreApplicationSettings();

   void DefineAtomGroup();
   void SelectAtomGroup();

   void TriggerLabelElementsGeneric();
public:
   void LoadApplicationSettings();
   // guess and set a size for the given sub-window in terms of the main window size.
   void SetDefaultDialogSize(QWidget *pWindow, double fDefaultVerticalScale = -1.);
protected:
   void dragEnterEvent(QDragEnterEvent *event); // override
   void dropEvent(QDropEvent *event); // override

   void AddPresetScript(QMenu *pMenu, QString FileName);
   void AssignScriptFileToAction(QAction *pAction, QString FileName);
   void closeEvent(QCloseEvent *event);

   // FNotifyEvents go here.
   void customEvent(QEvent *pEvent); // override
public:
   friend class IApplication;
   IApplication *ThisAsIApp();
//    void set_iso_surface_type_(QString const &IsoType, float fIsoValue);
   void load_files_(QStringList const &FileName);
   void close_files_();
   void notifyUser(FNotificationClass NotifyClass, QString const &Text, QString const &ExplicitTitle = QString());
protected:
   void MakeAtomGroupShortcuts();
   int iGetActionAtomGroupId();
};

// script interface
class IApplication : public FMainWindow {
   Q_OBJECT

   Q_PROPERTY(QStringList loaded_files READ get_loaded_files WRITE set_loaded_files)

   Q_PROPERTY(int orbital_color_scheme READ get_orbital_color_scheme WRITE set_orbital_color_scheme)
public slots:
   virtual void add_axes(QString Which, double AxisLength, QString Options);
   virtual void show_text(QString Caption, QString Text); // shows a ShowTextForm with a given text & caption
   virtual void notify(QString Text); // sets a text for the notification bar, as "notify" level.

   virtual void load_file(QString const &FileName); // = 0;
//    virtual void load_files(QString const &FileNames); // = 0;
   virtual void load_files(QScriptValue const &FileList);
   virtual void close_files();
   virtual void set_frame(int iFrame); // = 0;
   virtual int get_frame();
   virtual QString get_frame_name(); // = 0;
   virtual QStringList get_loaded_files();
   virtual void set_loaded_files(QStringList const &FileList); // calls close & load_files.
   // align molecules in space. 'Mode' denotes the atom weight; can
   // be "mass", "charge", "unity", or "iso-mass" (most common isotope mass instead of average isotope mass)
   virtual void orient_frames(QString const &Mode); // = 0;

   virtual void show_mo(int iMo, int cIsoPlus=-1, int cIsoMinus=-1); // = 0;
   virtual void hide_mo(int iMo); // = 0;
   virtual void hide_mos(); // = 0; // hide all MOs.
   virtual void rotate_mos_2x2(int iMo, int jMo, double Angle);

   virtual uint num_frames();
   virtual QObject* frame(int iFrame); // return frame #iFrame. return value is an IFrame* object.
   virtual QObject* frame(); // return current frame.

   // add/remove bond lines. Atom indices are 1-based here.
   // These here apply to all loaded frames.
   virtual void add_bond(int iAt, int jAt, QString const &Flags); // = 0;
   virtual void delete_bond(int iAt, int jAt); // = 0;
   virtual void reset_bonds(); // = 0; // reset all bonds to normal.

   virtual void set_atom_mode(int iAt, QString const &Mode); // = 0; // 0: hidden
   virtual void reset_atom_modes(); // = 0; // reset all atom states to normal.

   // exit the application at next opportunity.
   virtual void quit(); // = 0;

//    virtual FElementOptionsList &element_options();
//    virtual QObjectList element_options();
   virtual QObject* element_options(int iElem);
   // reset element options to defaults for that element.
   virtual void reset_element_options(int iElem);
   // note: this either returns or creates a override pPropertiesOverride object in the corresponding FAtomOptions object.
   virtual QObject* atom_options(int iAtom);
   // destroy property override object and make atom behave like other elements of its kind.
   virtual void reset_atom_options(int iAtom);
   virtual QObject* wf_options();

   virtual void update_views();

   virtual void define_atom_group(int iAtomGroup, QScriptValue const &AtomList);

   virtual int get_orbital_color_scheme();
   virtual void set_orbital_color_scheme(int);
//    virtual void set_iso_surface_type(QString const &IsoType, float fIsoValue); // = 0;

public: // here for technical reasons. Not part of script interface.
   IApplication( QWidget * parent = 0, Qt::WindowFlags flags = 0 )
      : FMainWindow(parent, flags)
   {}
   ~IApplication(); // does nothing---just to fix the vtable.
};


struct FAboutDialogImpl;

class FAboutDialog : public QDialog
{
    Q_OBJECT
public:
    Ui::AboutDialog *ui;
    FAboutDialogImpl *p;

    FAboutDialog(QWidget *parent = 0);
   ~FAboutDialog();
// protected:
//    void closeEvent(QCloseEvent *event);
};


extern int
   g_QtNotifyEventType; // goes to event::type (comes from QEvent::registerEventType, and will make something between QEvent::User and QEvent::MaxUser)
class FNotifyEvent : public QEvent
{
//     Q_OBJECT
   // ^-- QEvent doesn't inherit from Q_OBJECT, so doesn't have a Q_OBJECT initializer. Makes sense in retrospect...
public:
    FNotifyEvent(FNotificationClass NotifyClass, QString const &Message, QString const &Title = QString());
    virtual ~FNotifyEvent();
protected:
   FNotificationClass
      m_NotifyClass;
   QString
      m_Message;
   QString
      m_Title;
public:
   FNotificationClass notifyClass() const { return m_NotifyClass; }
   QString message() const { return m_Message; }
   QString title() const { return m_Title; }
};


// don't ask.
// class QFileDialog;
// class FSaveFileDirChangeProxy : public QObject
// {
//    Q_OBJECT
// public:
//    explicit FSaveFileDirChangeProxy(QFileDialog *pDialog);
//    ~FSaveFileDirChangeProxy();
// protected:
//    QFileDialog *m_pDialog;
// protected slots:
//    void selectFile(QString const &File);
// };

#endif // ORBVIEW_MAIN_H
