namespace ct {

// number of electrons in the orbitals of the respective shells,
// sorted by angular momentum.
float g_AtomicOccupations[1773] = {
   1,0,  // 1 Hydrogen: s
   2,0,  // 2 Helium: s
   2,1,0,  // 3 Lithium: s
   2,2,0,  // 4 Beryllium: s
   2,2,0, 1,0,  // 5 Boron: sp
   2,2,0, 2,0,  // 6 Carbon: sp
   2,2,0, 3,0,  // 7 Nitrogen: sp
   2,2,0, 4,0,  // 8 Oxygen: sp
   2,2,0, 5,0,  // 9 Fluorine: sp
   2,2,0, 6,0,  // 10 Neon: sp
   2,2,1,0, 6,0,  // 11 Sodium: sp
   2,2,2,0, 6,0,  // 12 Magnesium: sp
   2,2,2,0, 6,1,0,  // 13 Aluminium: sp
   2,2,2,0, 6,2,0,  // 14 Silicon: sp
   2,2,2,0, 6,3,0,  // 15 Phosphorus: sp
   2,2,2,0, 6,4,0,  // 16 Sulfur: sp
   2,2,2,0, 6,5,0,  // 17 Chlorine: sp
   2,2,2,0, 6,6,0,  // 18 Argon: sp
   2,2,2,1,0, 6,6,0,  // 19 Potassium: sp
   2,2,2,2,0, 6,6,0,  // 20 Calcium: sp
   2,2,2,2,0, 6,6,0, 1,0,  // 21 Scandium: spd
   2,2,2,2,0, 6,6,0, 2,0,  // 22 Titanium: spd
   2,2,2,2,0, 6,6,0, 3,0,  // 23 Vanadium: spd
   2,2,2,1,0, 6,6,0, 5,0,  // 24 Chromium: spd
   2,2,2,2.00,0, 6,6,0, 5.00,0,  // 25 Manganese: spd
   2,2,2,2,0, 6,6,0, 6,0,  // 26 Iron: spd
   2,2,2,2,0, 6,6,0, 7,0,  // 27 Cobalt: spd
   2,2,2,2,0, 6,6,0, 8,0,  // 28 Nickel: spd
   2,2,2,1,0, 6,6,0, 10,0,  // 29 Copper: spd
   2,2,2,2,0, 6,6,0, 10,0,  // 30 Zinc: spd
   2,2,2,2,0, 6,6,1,0, 10,0,  // 31 Gallium: spd
   2,2,2,2,0, 6,6,2,0, 10,0,  // 32 Germanium: spd
   2,2,2,2,0, 6,6,3,0, 10,0,  // 33 Arsenic: spd
   2,2,2,2,0, 6,6,4,0, 10,0,  // 34 Selenium: spd
   2,2,2,2,0, 6,6,5,0, 10,0,  // 35 Bromine: spd
   2,2,2,2,0, 6,6,6,0, 10,0,  // 36 Krypton: spd
   2,2,2,2,1,0, 6,6,6,0, 10,0,  // 37 Rubidium: spd
   2,2,2,2,2,0, 6,6,6,0, 10,0,  // 38 Strontium: spd
   2,2,2,2,2,0, 6,6,6,0, 10,1,0,  // 39 Yttrium: spd
   2,2,2,2,2,0, 6,6,6,0, 10,2,0,  // 40 Zirconium: spd
   2,2,2,2,1,0, 6,6,6,0, 10,4,0,  // 41 Niobium: spd
   2,2,2,2,1,0, 6,6,6,0, 10,5,0,  // 42 Molybdenum: spd
   2,2,2,2,2,0, 6,6,6,0, 10,5,0,  // 43 Technetium: spd
   2,2,2,2,1,0, 6,6,6,0, 10,7,0,  // 44 Ruthenium: spd
   2,2,2,2,1,0, 6,6,6,0, 10,8,0,  // 45 Rhodium: spd
//    2,2,2,2,0,0, 6,6,6,0, 10,10,0,  // 46 Palladium: spd
// ^-- UPDATE (cgk, 2022-01-22): gna... this gets read as:
//
//     Neutral atom shell occupation:  s[2222] p[] d[666] f[aa]
//
//     Thanks to Prof Russel Hughes from Darthmouth for providing a simple case
//     file which shows the problem. Apparently I indeed never once did a single
//     calculation involving an SCF guess which involved Palladium.
//
//     Now behold my super elegant fix to this problem:
   2,2,2,2,0.000001,0, 6,6,5.999999,0, 10,10,0,  // 46 Palladium: spd
// ^-- why not just "2,2,2,2,0, 6,6,6,0, 10,10,0"? Because that would
//     require reindexing the table (g_AtomicOccupationRanges below),
//     and I do not think I've touched the script which made it for
//     a decade. I should still have it *somewhere*, I guess...
   2,2,2,2,1,0, 6,6,6,0, 10,10,0,  // 47 Silver: spd
   2,2,2,2,2,0, 6,6,6,0, 10,10,0,  // 48 Cadmium: spd
   2,2,2,2,2,0, 6,6,6,1,0, 10,10,0,  // 49 Indium: spd
   2,2,2,2,2,0, 6,6,6,2,0, 10,10,0,  // 50 Tin: spd
   2,2,2,2,2,0, 6,6,6,3,0, 10,10,0,  // 51 Antimony: spd
   2,2,2,2,2,0, 6,6,6,4,0, 10,10,0,  // 52 Tellurium: spd
   2,2,2,2,2,0, 6,6,6,5,0, 10,10,0,  // 53 Iodine: spd
   2,2,2,2,2,0, 6,6,6,6,0, 10,10,0,  // 54 Xenon: spd
   2,2,2,2,2,1,0, 6,6,6,6,0, 10,10,0,  // 55 Caesium: spd
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0,  // 56 Barium: spd
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,1,0,  // 57 Lanthanum: spd
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,1,0, 1,0,  // 58 Cerium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 3,0,  // 59 Praseodymium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 4,0,  // 60 Neodymium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 5,0,  // 61 Promethium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 6,0,  // 62 Samarium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 7,0,  // 63 Europium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,1,0, 7,0,  // 64 Gadolinium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 9,0,  // 65 Terbium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 10,0,  // 66 Dysprosium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 11,0,  // 67 Holmium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 12,0,  // 68 Erbium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 13,0,  // 69 Thulium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,0, 14,0,  // 70 Ytterbium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,1,0, 14,0,  // 71 Lutetium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,2,0, 14,0,  // 72 Hafnium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,3,0, 14,0,  // 73 Tantalum: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,4,0, 14,0,  // 74 Tungsten: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,5,0, 14,0,  // 75 Rhenium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,6,0, 14,0,  // 76 Osmium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,7,0, 14,0,  // 77 Iridium: spdf
   2,2,2,2,2,1,0, 6,6,6,6,0, 10,10,9,0, 14,0,  // 78 Platinum: spdf
   2,2,2,2,2,1,0, 6,6,6,6,0, 10,10,10,0, 14,0,  // 79 Gold: spdf
   2,2,2,2,2,2,0, 6,6,6,6,0, 10,10,10,0, 14,0,  // 80 Mercury: spdf
   2,2,2,2,2,2,0, 6,6,6,6,1,0, 10,10,10,0, 14,0,  // 81 Thallium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,2,0, 10,10,10,0, 14,0,  // 82 Lead: spdf
   2,2,2,2,2,2,0, 6,6,6,6,3,0, 10,10,10,0, 14,0,  // 83 Bismuth: spdf
   2,2,2,2,2,2,0, 6,6,6,6,4,0, 10,10,10,0, 14,0,  // 84 Polonium: spdf
   2,2,2,2,2,2,0, 6,6,6,6,5,0, 10,10,10,0, 14,0,  // 85 Astatine: spdf
   2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,0,  // 86 Radon: spdf
   2,2,2,2,2,2,1,0, 6,6,6,6,6,0, 10,10,10,0, 14,0,  // 87 Francium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,0,  // 88 Radium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,1,0, 14,0,  // 89 Actinium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,2,0, 14,0,  // 90 Thorium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,1,0, 14,2,0,  // 91 Protactinium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,1,0, 14,3,0,  // 92 Uranium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,1,0, 14,4,0,  // 93 Neptunium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,6,0,  // 94 Plutonium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,7,0,  // 95 Americium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,1,0, 14,7,0,  // 96 Curium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,9,0,  // 97 Berkelium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,10,0,  // 98 Californium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,11,0,  // 99 Einsteinium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,12,0,  // 100 Fermium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,13,0,  // 101 Mendelevium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,0, 14,14,0,  // 102 Nobelium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,1,0, 14,14,0,  // 103 Lawrencium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,2,0, 14,14,0,  // 104 Rutherfordium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,3,0, 14,14,0,  // 105 Dubnium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,4,0, 14,14,0,  // 106 Seaborgium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,5,0, 14,14,0,  // 107 Bohrium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,6,0, 14,14,0,  // 108 Hassium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,7,0, 14,14,0,  // 109 Meitnerium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,8,0, 14,14,0,  // 110 Darmstadtium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,9,0, 14,14,0,  // 111 Roentgenium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,0, 10,10,10,10,0, 14,14,0,  // 112 Copernicium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,1,0, 10,10,10,10,0, 14,14,0,  // 113 Ununtrium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,2,0, 10,10,10,10,0, 14,14,0,  // 114 Ununquadium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,3,0, 10,10,10,10,0, 14,14,0,  // 115 Ununpentium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,4,0, 10,10,10,10,0, 14,14,0,  // 116 Ununhexium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,5,0, 10,10,10,10,0, 14,14,0,  // 117 Ununseptium: spdf
   2,2,2,2,2,2,2,0, 6,6,6,6,6,6,0, 10,10,10,10,0, 14,14,0,  // 118 Ununoctium: spd
};

// index into g_AtomicOccupations.
short g_AtomicOccupationRanges[118][2] = {
   {0, 2}, // Hydrogen
   {2, 4}, // Helium
   {4, 7}, // Lithium
   {7, 10}, // Beryllium
   {10, 15}, // Boron
   {15, 20}, // Carbon
   {20, 25}, // Nitrogen
   {25, 30}, // Oxygen
   {30, 35}, // Fluorine
   {35, 40}, // Neon
   {40, 46}, // Sodium
   {46, 52}, // Magnesium
   {52, 59}, // Aluminium
   {59, 66}, // Silicon
   {66, 73}, // Phosphorus
   {73, 80}, // Sulfur
   {80, 87}, // Chlorine
   {87, 94}, // Argon
   {94, 102}, // Potassium
   {102, 110}, // Calcium
   {110, 120}, // Scandium
   {120, 130}, // Titanium
   {130, 140}, // Vanadium
   {140, 150}, // Chromium
   {150, 160}, // Manganese
   {160, 170}, // Iron
   {170, 180}, // Cobalt
   {180, 190}, // Nickel
   {190, 200}, // Copper
   {200, 210}, // Zinc
   {210, 221}, // Gallium
   {221, 232}, // Germanium
   {232, 243}, // Arsenic
   {243, 254}, // Selenium
   {254, 265}, // Bromine
   {265, 276}, // Krypton
   {276, 288}, // Rubidium
   {288, 300}, // Strontium
   {300, 313}, // Yttrium
   {313, 326}, // Zirconium
   {326, 339}, // Niobium
   {339, 352}, // Molybdenum
   {352, 365}, // Technetium
   {365, 378}, // Ruthenium
   {378, 391}, // Rhodium
   {391, 404}, // Palladium
   {404, 417}, // Silver
   {417, 430}, // Cadmium
   {430, 444}, // Indium
   {444, 458}, // Tin
   {458, 472}, // Antimony
   {472, 486}, // Tellurium
   {486, 500}, // Iodine
   {500, 514}, // Xenon
   {514, 529}, // Caesium
   {529, 544}, // Barium
   {544, 560}, // Lanthanum
   {560, 578}, // Cerium
   {578, 595}, // Praseodymium
   {595, 612}, // Neodymium
   {612, 629}, // Promethium
   {629, 646}, // Samarium
   {646, 663}, // Europium
   {663, 681}, // Gadolinium
   {681, 698}, // Terbium
   {698, 715}, // Dysprosium
   {715, 732}, // Holmium
   {732, 749}, // Erbium
   {749, 766}, // Thulium
   {766, 783}, // Ytterbium
   {783, 801}, // Lutetium
   {801, 819}, // Hafnium
   {819, 837}, // Tantalum
   {837, 855}, // Tungsten
   {855, 873}, // Rhenium
   {873, 891}, // Osmium
   {891, 909}, // Iridium
   {909, 927}, // Platinum
   {927, 945}, // Gold
   {945, 963}, // Mercury
   {963, 982}, // Thallium
   {982, 1001}, // Lead
   {1001, 1020}, // Bismuth
   {1020, 1039}, // Polonium
   {1039, 1058}, // Astatine
   {1058, 1077}, // Radon
   {1077, 1097}, // Francium
   {1097, 1117}, // Radium
   {1117, 1138}, // Actinium
   {1138, 1159}, // Thorium
   {1159, 1181}, // Protactinium
   {1181, 1203}, // Uranium
   {1203, 1225}, // Neptunium
   {1225, 1246}, // Plutonium
   {1246, 1267}, // Americium
   {1267, 1289}, // Curium
   {1289, 1310}, // Berkelium
   {1310, 1331}, // Californium
   {1331, 1352}, // Einsteinium
   {1352, 1373}, // Fermium
   {1373, 1394}, // Mendelevium
   {1394, 1415}, // Nobelium
   {1415, 1437}, // Lawrencium
   {1437, 1459}, // Rutherfordium
   {1459, 1481}, // Dubnium
   {1481, 1503}, // Seaborgium
   {1503, 1525}, // Bohrium
   {1525, 1547}, // Hassium
   {1547, 1569}, // Meitnerium
   {1569, 1591}, // Darmstadtium
   {1591, 1613}, // Roentgenium
   {1613, 1635}, // Copernicium
   {1635, 1658}, // Ununtrium
   {1658, 1681}, // Ununquadium
   {1681, 1704}, // Ununpentium
   {1704, 1727}, // Ununhexium
   {1727, 1750}, // Ununseptium
   {1750, 1773}  // Ununoctium
};

} // namespace ct
