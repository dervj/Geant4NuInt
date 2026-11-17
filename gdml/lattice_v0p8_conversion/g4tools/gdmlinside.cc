// gdmlinside.cc -- overlap/gap check by point sampling, independent of
// G4VSolid::GetPointOnSurface (which fails for the deep pyg4ometry booleans).
// Random points in a box; for each point, at every level of the hierarchy,
// count the sibling daughters whose solid has the point kInside.
//   >1 => overlap (reported per pair);  0 inside a cell => gap in the cell.
//   gdmlinside file.gdml N x0 x1 y0 y1 z0 z1 (mm) [seed]
#include <G4AffineTransform.hh>
#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <string>
int main(int argc, char** argv)
{
    G4GDMLParser parser;
    parser.Read(argv[1], false);
    long N = std::atol(argv[2]);
    double b[6];
    for (int i = 0; i < 6; ++i) b[i] = std::atof(argv[3 + i]);
    std::mt19937_64 rng(argc > 9 ? std::atol(argv[9]) : 1);
    std::uniform_real_distribution<double> U(0, 1);
    G4LogicalVolume* world = parser.GetWorldVolume()->GetLogicalVolume();
    std::map<std::string, long> pairs, gaps;
    long nworldgap = 0;
    for (long n = 0; n < N; ++n)
    {
        G4ThreeVector p(b[0] + (b[1] - b[0]) * U(rng), b[2] + (b[3] - b[2]) * U(rng),
                        b[4] + (b[5] - b[4]) * U(rng));
        G4LogicalVolume* mother = world;
        G4ThreeVector lp = p;
        for (int depth = 0; depth < 10; ++depth)
        {
            G4VPhysicalVolume* hit = nullptr;
            G4ThreeVector hitp;
            for (std::size_t i = 0; i < mother->GetNoDaughters(); ++i)
            {
                G4VPhysicalVolume* d = mother->GetDaughter(i);
                G4AffineTransform t(d->GetRotation(), d->GetTranslation());
                G4ThreeVector q = t.Inverse().TransformPoint(lp);
                if (d->GetLogicalVolume()->GetSolid()->Inside(q) == kInside)
                {
                    if (hit) {
                        std::string a = hit->GetName(), c = d->GetName();
                        if (c < a) std::swap(a, c);
                        if (pairs[a + " | " + c]++ == 0)
                            std::printf("OVL %s | %s in %s at (%.1f,%.1f,%.1f)\n", a.c_str(), c.c_str(),
                                        mother->GetName().c_str(), p.x(), p.y(), p.z());
                    } else { hit = d; hitp = q; }
                }
            }
            if (!hit) {
                if (mother != world && mother->GetNoDaughters() > 0) {
                    if (gaps[mother->GetName()]++ < 3)
                        std::printf("GAP in %s at (%.2f,%.2f,%.2f)\n", mother->GetName().c_str(), p.x(), p.y(), p.z());
                } else if (mother == world) ++nworldgap;
                break;
            }
            mother = hit->GetLogicalVolume();
            lp = hitp;
        }
    }
    std::printf("SUMMARY N=%ld overlapping_pairs=%zu cells_with_gaps=%zu world_gap_points=%ld\n", N, pairs.size(),
                gaps.size(), nworldgap);
    for (auto& kv : pairs) std::printf("  pair %s : %ld pts\n", kv.first.c_str(), kv.second);
    for (auto& kv : gaps) std::printf("  gap %s : %ld pts\n", kv.first.c_str(), kv.second);
}
