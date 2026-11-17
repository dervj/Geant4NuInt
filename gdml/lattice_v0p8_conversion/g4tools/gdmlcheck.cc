// gdmlcheck.cc -- Geant4-only sanity check of a GDML file before handing it
// to Celeritas: world extent, each world daughter's bounding box (in world
// coordinates) vs the world box, and optionally Geant4's overlap check.
//   gdmlcheck file.gdml [overlap_points=0] [only_name_substring]
#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4Material.hh>
#include <G4AffineTransform.hh>
#include <G4ios.hh>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: %s file.gdml [npts] [name]\n", argv[0]); return 2; }
    int npts = argc > 2 ? std::atoi(argv[2]) : 0;
    std::string only = argc > 3 ? argv[3] : "";
    G4GDMLParser parser;
    parser.SetOverlapCheck(false);
    parser.Read(argv[1], false);
    G4VPhysicalVolume* world = parser.GetWorldVolume();
    G4LogicalVolume* wlv = world->GetLogicalVolume();
    G4ThreeVector wlo, whi;
    wlv->GetSolid()->BoundingLimits(wlo, whi);
    std::printf("WORLD %s solid=%s type=%s lo=(%.1f,%.1f,%.1f) hi=(%.1f,%.1f,%.1f) mm, %zu daughters\n",
                world->GetName().c_str(), wlv->GetSolid()->GetName().c_str(),
                wlv->GetSolid()->GetEntityType().c_str(), wlo.x(), wlo.y(), wlo.z(),
                whi.x(), whi.y(), whi.z(), wlv->GetNoDaughters());
    int n_out = 0, n_partial = 0, n_ovl = 0;
    for (std::size_t i = 0; i < wlv->GetNoDaughters(); ++i)
    {
        G4VPhysicalVolume* pv = wlv->GetDaughter(i);
        G4VSolid* s = pv->GetLogicalVolume()->GetSolid();
        G4ThreeVector lo, hi;
        s->BoundingLimits(lo, hi);
        G4AffineTransform t(pv->GetRotation(), pv->GetTranslation());
        G4ThreeVector glo(1e30, 1e30, 1e30), ghi(-1e30, -1e30, -1e30);
        for (int c = 0; c < 8; ++c)
        {
            G4ThreeVector p((c & 1) ? hi.x() : lo.x(), (c & 2) ? hi.y() : lo.y(), (c & 4) ? hi.z() : lo.z());
            // G4AffineTransform(rot, tlate) built from a placement maps local->mother
            // with the inverse rotation convention of G4PVPlacement
            G4ThreeVector g = (pv->GetRotation() ? pv->GetRotation()->inverse() * p : p) + pv->GetTranslation();
            glo = G4ThreeVector(std::min(glo.x(), g.x()), std::min(glo.y(), g.y()), std::min(glo.z(), g.z()));
            ghi = G4ThreeVector(std::max(ghi.x(), g.x()), std::max(ghi.y(), g.y()), std::max(ghi.z(), g.z()));
        }
        bool disjoint = ghi.x() < wlo.x() || glo.x() > whi.x() || ghi.y() < wlo.y() || glo.y() > whi.y()
                        || ghi.z() < wlo.z() || glo.z() > whi.z();
        bool inside = glo.x() >= wlo.x() && ghi.x() <= whi.x() && glo.y() >= wlo.y() && ghi.y() <= whi.y()
                      && glo.z() >= wlo.z() && ghi.z() <= whi.z();
        if (disjoint) ++n_out; else if (!inside) ++n_partial;
        bool ovl = false;
        if (npts > 0 && (only.empty() || pv->GetName().find(only) != std::string::npos))
        {
            ovl = pv->CheckOverlaps(npts, 0., false, 1);
            if (ovl) ++n_ovl;
        }
        std::printf("PV %-28s %-12s %-10s bbox lo=(%.0f,%.0f,%.0f) hi=(%.0f,%.0f,%.0f) %s%s\n",
                    pv->GetName().c_str(), pv->GetLogicalVolume()->GetMaterial()->GetName().c_str(),
                    s->GetEntityType().c_str(), glo.x(), glo.y(), glo.z(), ghi.x(), ghi.y(), ghi.z(),
                    disjoint ? "OUTSIDE_WORLD" : (inside ? "in" : "PROTRUDES"),
                    npts > 0 ? (ovl ? " OVERLAP" : " ovl-ok") : "");
    }
    std::printf("SUMMARY daughters=%zu fully_outside_world=%d protruding=%d overlaps=%d (npts=%d)\n",
                wlv->GetNoDaughters(), n_out, n_partial, n_ovl, npts);
    return 0;
}
