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

#ifndef SCRIPT_INTERFACE_H
#define SCRIPT_INTERFACE_H

#include "Iv.h"
#include "IvMain.h"
#include "IvView3D.h"
#include "IvDocument.h"
#include <QString>

void ExecScript(IApplication *app, IView3d *view, QString const &ScriptText, QString const &ScriptName);
void ExecScript(IApplication *app, IView3d *view, QString const &FileName);
QString LoadTextFileViaQt(QString const &FileName);

#endif // SCRIPT_INTERFACE_H
