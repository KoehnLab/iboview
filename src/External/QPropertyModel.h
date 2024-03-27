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

// QPropertyModel
// - a class for easily turning any QObject-derived subclass with properties into a one-row model
//
// Copyright 2013 - Harvey Chapman <hchapman@3gfp.com>
// Source: https://gist.github.com/sr105/7955969
// License:
//   This work is licensed under the Creative Commons Attribution-ShareAlike
//   4.0 International License. To view a copy of this license, visit
//   http://creativecommons.org/licenses/by-sa/4.0/deed.en_US.
//
// It's not required, but I'd appreciate it if any improvements were e-mailed
// back to me so I can share them with others. This code is specifically not
// GPL-like so you can use it commercially without worrying about it tainting
// the rest of your proprietary code.
// -- Harvey

// Notes by cgk:
// - Taken from https://gist.github.com/sr105/7955969
// - See also: http://stackoverflow.com/questions/18793735/connect-a-signal-to-the-slot-of-a-qmetaproperty
//   (which is what I originally wanted to do, but failed to find a viable way, just as the
//    author of that question)
// - Various slight changes (in particular, un-staticify the property name list/map)

#ifndef QPROPERTYMODEL_H
#define QPROPERTYMODEL_H

#include <QAbstractItemModel>
#include <QStringList>
#include <QDataWidgetMapper>
#include <QMap>
#include <QMetaProperty>
#include <QAction>

void LinkPropertyWidgets(QObject *pTarget, QWidget *pWidgetContainer, char const *pPropertyKeyName);

class QPropertyModel;

// Convenience class that exposes the public methods of QPropertyModel
// without requiring casting.
class QPropertyDataWidgetMapper : public QDataWidgetMapper
{
    Q_OBJECT
public:
    QPropertyDataWidgetMapper(QObject *parent = 0) : QDataWidgetMapper(parent) {}

    // QDataWidgetMapper::model() re-written to return QPropertyDataWidgetMapper
    QPropertyModel *model() const;
    // For convenience, these automatically convert "property" into column numbers
    void addMapping(QWidget *widget, QString property);
    void addMapping(QWidget *widget, QString property, const QByteArray &propertyName);
//     void addMapping(QAction *action, QString property);
    // Pass-thru methods to QDataWidgetMapper
    void addMapping(QWidget *widget, int section);
    void addMapping(QWidget *widget, int section, const QByteArray &propertyName);
//     void addMapping(QAction *action, int section);
};

// QPropertyModel creates a single row data model consisting of columns mapping
// to properties in a QObject. The column list can be retrieved as a QStringList,
// and a method exists to convert the property names to column numbers.
class QPropertyModel : public QAbstractItemModel
{
    Q_OBJECT
public:
    explicit QPropertyModel(QObject *source, QObject *parent = 0);
    ~QPropertyModel();

    // Return a QPropertyDataWidgetMapper wrapping a new instance of this class.
    static QPropertyDataWidgetMapper *newMapper(QObject *source, QObject *parent = 0);
    // Return a QPropertyDataWidgetMapper wrapping this existing instance
    QPropertyDataWidgetMapper *newMapper();

    QStringList propertyNames() const;
    int columnForProperty(QString name) const;
    QMap<int, QMetaProperty> properties() const;

protected:
    void connectToPropertyNotifySignals();

    // Pointer to our data source
    QObject *_source;

    mutable QStringList _names;
    mutable QMap<int, QMetaProperty> _properties;
protected slots:
    void columnChanged(int column);

    // Required virtual function implementations. They mostly map
    // directly to the (getItem/setItem/itemChanged) methods above.
public:
    // read & write data
    virtual QVariant data(const QModelIndex &index, int role) const;
    virtual bool setData(const QModelIndex &index, const QVariant &value, int role);

    // returns the number of properties in _source
    virtual int columnCount(const QModelIndex &parent) const;

    // all hard-coded simple implementations
    virtual int rowCount(const QModelIndex &parent = QModelIndex()) const;
    virtual Qt::ItemFlags flags(const QModelIndex &index) const;
    virtual QModelIndex parent(const QModelIndex &child) const;
    virtual QModelIndex index(int row, int column, const QModelIndex &) const;

    // Helper method to make virtual methods easier to code
    virtual bool hasIndex(const QModelIndex &index) const;
};


// Until we can come up with something more clever, this little class allows
// us to connect each signal in a single QObject to a single slot using
// QSignalMapper to pass information to us about which signal was sent.
// QSignalMapper maps Objects to data. All of our signals come from the same
// object, so that won't work. However, if we create a SignalObject as a
// forwareder for each signal, now we have a unique object for each signal
// that QSignalMapper can work with.
class SignalForwarder: public QObject
{
    Q_OBJECT
public:
    SignalForwarder(QObject *parent = 0) : QObject(parent) {}
signals:
    void forward();
};

#endif // QPROPERTYMODEL_H
