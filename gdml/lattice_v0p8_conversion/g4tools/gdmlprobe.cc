// gdmlprobe.cc -- locate many points in a GDML world.
//   gdmlprobe file.gdml points.txt   (points.txt: "x y z" in mm per line)
// prints: idx  LV  material  motherLV  depth
#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Navigator.hh>
#include <G4TouchableHistory.hh>
#include <G4VPhysicalVolume.hh>
#include <cstdio>
#include <fstream>
int main(int argc, char** argv)
{
    G4GDMLParser parser;
    parser.Read(argv[1], false);
    G4Navigator nav;
    nav.SetWorldVolume(parser.GetWorldVolume());
    std::ifstream in(argv[2]);
    double x, y, z;
    long i = 0;
    while (in >> x >> y >> z)
    {
        nav.LocateGlobalPointAndSetup({x, y, z}, nullptr, false, true);
        G4TouchableHistory* t = nav.CreateTouchableHistory();
        int d = t->GetHistoryDepth();
        G4VPhysicalVolume* pv = t->GetVolume(0);
        const char* lv = pv ? pv->GetLogicalVolume()->GetName().c_str() : "<outside>";
        const char* mat = pv ? pv->GetLogicalVolume()->GetMaterial()->GetName().c_str() : "-";
        const char* mo = (pv && d >= 1) ? t->GetVolume(1)->GetLogicalVolume()->GetName().c_str() : "-";
        std::printf("%ld %s %s %s %d\n", i++, lv, mat, mo, d);
        delete t;
    }
}
