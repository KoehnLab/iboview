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

#include "QPropertyModel.h"

#include <QDataWidgetMapper>
#include <QMetaProperty>
#include <QSignalMapper>
#include <QDebug>

QPropertyModel::QPropertyModel(QObject *source, QObject *parent) :
    QAbstractItemModel(parent),
    _source(source)
{
    connectToPropertyNotifySignals();
}

QPropertyModel::~QPropertyModel()
{
}

QPropertyDataWidgetMapper *QPropertyModel::newMapper(QObject *source, QObject *parent)
{
    QPropertyDataWidgetMapper *_mapper = new QPropertyDataWidgetMapper(parent);
    _mapper->setModel(new QPropertyModel(source, parent));
    return _mapper;
}

QPropertyDataWidgetMapper *QPropertyModel::newMapper()
{
    QPropertyDataWidgetMapper *_mapper = new QPropertyDataWidgetMapper(this);
    _mapper->setModel(this);
    return _mapper;
}

QStringList QPropertyModel::propertyNames() const
{
//     static QStringList names;
// cgk: ^- hmpf. Why would one make those static? Took >1hr to find this...
    if (_names.size())
        return _names;

    foreach (const QMetaProperty &p, properties().values())
        _names << QString::fromLatin1(p.name());
    return _names;
}

int QPropertyModel::columnForProperty(QString name) const
{
    QMap<int, QMetaProperty> props = properties();
    foreach (int index, props.keys())
        if (props[index].name() == name)
            return index;
    qDebug("No property \"%s\" found!", qPrintable(name));
    return -1;
}

QMap<int, QMetaProperty> QPropertyModel::properties() const
{
    // TODO: should we filter out properties with no NOTIFY? Otherwise, how will
    //       we know when they change? Maybe it doesn't matter. Just let the user
    //       any QObject they want and deal with the consequences of the features
    //       that their QObject properties support.
//     static QMap<int, QMetaProperty> properties;
    if (_properties.size())
        return _properties;

    // Save a map of properties
    const QMetaObject* metaObject = _source->metaObject();
    // Start at 0 to get all inherited properties, too. Start at the offset for just this
    // subclass.
    // for(int i = metaObject->propertyOffset(); i < metaObject->propertyCount(); ++i)
    for(int i = 0; i < metaObject->propertyCount(); ++i)
         _properties.insert(i, metaObject->property(i));
    return _properties;
}

void QPropertyModel::connectToPropertyNotifySignals()
{
    QMap<int, QMetaProperty> props = properties();
    QSignalMapper *mapper = new QSignalMapper(this);
    foreach (int index, props.keys()) {
        if (!props[index].hasNotifySignal())
            continue;
        // It's difficult to map all signals from a single object to
        // a single slot with an identifiable piece of information.
        // Ideas:
        //   - http://stackoverflow.com/questions/10805174/qobject-generic-signal-handler
        //   - a dummy SignalForwarder class instance for each signal for QSignalMapper.
        SignalForwarder *sf = new SignalForwarder(this);
#if QT_VERSION <= 0x050000
        connect(_source, QByteArray("2") + props[index].notifySignal().signature(), sf, SIGNAL(forward()));
#else
        connect(_source, QByteArray("2") + props[index].notifySignal().methodSignature(), sf, SIGNAL(forward()));
#endif
        connect(sf, SIGNAL(forward()), mapper, SLOT(map()));
        mapper->setMapping(sf, index);
    }
    connect(mapper, SIGNAL(mapped(int)), this, SLOT(columnChanged(int)));
}

void QPropertyModel::columnChanged(int column)
{
    //qDebug("columnChanged(%d)", column);
    QModelIndex index = createIndex(0, column);
    emit dataChanged(index, index);
}

QVariant QPropertyModel::data(const QModelIndex& index, int role) const {
    //qDebug() << "data(" << index << "," << role << ")";
    if (!hasIndex(index))
        return QVariant();
    if (role == Qt::DisplayRole || role == Qt::EditRole) // <- cgk: added role checks here.
        return _source->property(properties().value(index.column()).name());
    if (role == Qt::ToolTipRole)
        return QString("WHEEEEE!!!!!");
    return QVariant();
}

bool QPropertyModel::setData(const QModelIndex &index, const QVariant &value, int role) {
    //qDebug() << "setData(" << index << "," << value << "," << role << ")";
    if (!hasIndex(index) || role != Qt::EditRole)
        return QAbstractItemModel::setData(index, value, role);

    QMetaProperty mp = properties().value(index.column());
    bool rc = _source->setProperty(mp.name(), value);
    if (rc && !mp.hasNotifySignal())
        // property doesn't support NOTIFY, emit dataChanged()
        emit dataChanged(index, index);
    return rc;
}

int QPropertyModel::columnCount(const QModelIndex &/*parent*/) const {
    return properties().size();
}

int QPropertyModel::rowCount(const QModelIndex &/*parent*/) const {
    return 1;
}

Qt::ItemFlags QPropertyModel::flags(const QModelIndex &index) const {
    QMetaProperty mp = properties().value(index.column());
    return QAbstractItemModel::flags(index)
            | Qt::ItemIsSelectable
            | Qt::ItemIsEnabled
            | (mp.isWritable() ? Qt::ItemIsEditable : Qt::NoItemFlags);
}

QModelIndex QPropertyModel::parent(const QModelIndex &/*child*/) const {
    return QModelIndex();
}

QModelIndex QPropertyModel::index(int row, int column, const QModelIndex &parent) const {
    if (QAbstractItemModel::hasIndex(row, column, parent))
        return createIndex(row, column);
    return QModelIndex();
}

bool QPropertyModel::hasIndex(const QModelIndex &index) const {
    return QAbstractItemModel::hasIndex(index.row(), index.column(), index.parent());
}

QPropertyModel *QPropertyDataWidgetMapper::model() const
{
    return qobject_cast<QPropertyModel*>(QDataWidgetMapper::model());
}

void QPropertyDataWidgetMapper::addMapping(QWidget *widget, QString property)
{
    if (model() && model()->columnForProperty(property) >= 0)
        QDataWidgetMapper::addMapping(widget, model()->columnForProperty(property));
}

void QPropertyDataWidgetMapper::addMapping(QWidget *widget, QString property, const QByteArray &propertyName)
{
    if (model() && model()->columnForProperty(property) >= 0)
        QDataWidgetMapper::addMapping(widget, model()->columnForProperty(property), propertyName);
}

void QPropertyDataWidgetMapper::addMapping(QWidget *widget, int section)
{
    QDataWidgetMapper::addMapping(widget, section);
}

void QPropertyDataWidgetMapper::addMapping(QWidget *widget, int section, const QByteArray &propertyName)
{
    QDataWidgetMapper::addMapping(widget, section, propertyName);
}

// void QPropertyDataWidgetMapper::addMapping(QAction *action, QString property)
// {
//     if (model() && model()->columnForProperty(property) >= 0)
//         QDataWidgetMapper::addMapping(action, model()->columnForProperty(property));
// }
//
//
// void QPropertyDataWidgetMapper::addMapping(QAction *action, int section)
// {
//     QDataWidgetMapper::addMapping(action, section);
// }
