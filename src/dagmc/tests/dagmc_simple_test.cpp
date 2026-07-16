#include <gtest/gtest.h>

#include <iostream>
#include <limits>
#include <vector>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"

using namespace moab;

using moab::DagMC;

moab::DagMC* DAG;

static const char input_file[] = "test_geom.h5m";

static Range geometry_vertices(DagMC& dagmc) {
  Range surfaces, triangles, vertices;
  EXPECT_EQ(dagmc.geom_tool()->get_gsets_by_dimension(2, surfaces), MB_SUCCESS);
  for (EntityHandle surface : surfaces)
    EXPECT_EQ(dagmc.moab_instance()->get_entities_by_type(surface, MBTRI,
                                                          triangles, true),
              MB_SUCCESS);
  EXPECT_EQ(dagmc.moab_instance()->get_adjacencies(triangles, 0, false,
                                                   vertices, Interface::UNION),
            MB_SUCCESS);
  EXPECT_LT(vertices.size(), 3 * triangles.size());
  return vertices;
}

static std::vector<double> geometry_coordinates(DagMC& dagmc) {
  Range vertices = geometry_vertices(dagmc);
  std::vector<double> coordinates(3 * vertices.size());
  EXPECT_EQ(dagmc.moab_instance()->get_coords(vertices, coordinates.data()),
            MB_SUCCESS);
  return coordinates;
}

static size_t count_obb_roots(Interface* mbi) {
  Tag tag;
  if (mbi->tag_get_handle("OBB_ROOT", tag) != MB_SUCCESS) return 0;
  Range tagged_sets;
  EXPECT_EQ(mbi->get_entities_by_type_and_tag(0, MBENTITYSET, &tag, NULL, 1,
                                              tagged_sets),
            MB_SUCCESS);
  return tagged_sets.size();
}

static double ray_distance(DagMC& dagmc) {
  double origin[3] = {0.0, 0.0, 0.0};
  double direction[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surface;
  double distance = 0.0;
  EXPECT_EQ(dagmc.ray_fire(dagmc.entity_by_index(3, 1), origin, direction,
                           next_surface, distance),
            MB_SUCCESS);
  return distance;
}

class DagmcSimpleTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(DagmcSimpleTest, dagmc_load_file) {
  DAG = new DagMC();
  ErrorCode rval = DAG->load_file(input_file);  // open the Dag file
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc) {
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  // make new dagmc into that moab
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);

  ErrorCode rval;

  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_via_moab) {
  /* 2 - Test with external moab, load file in MOAB*/
  // load the file into moab rather than dagmc
  ErrorCode rval;

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  rval = mbi->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);
  rval = dagmc->load_existing_contents();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_internal) {
  /* 3 - Test with internal moab, load file in DAG*/
  // make new dagmc into that moab
  ErrorCode rval;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_build_obb) {
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  ErrorCode rval;

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  // make new dagmc into that moab
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);

  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_via_moab_build_obb) {
  /* 2 - Test with external moab, load file in MOAB*/
  // load the file into moab rather than dagmc
  ErrorCode rval;

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  rval = mbi->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);
  rval = dagmc->load_existing_contents();
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_internal_build_obb) {
  /* 3 - Test with internal moab, load file in DAG*/
  // make new dagmc into that moab
  ErrorCode rval;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_length_multiplier) {
  DagMC baseline;
  ASSERT_EQ(baseline.load_file(input_file), MB_SUCCESS);
  EXPECT_THROW(baseline.set_length_multiplier(2.0), std::logic_error);
  std::vector<double> baseline_coords = geometry_coordinates(baseline);
  ASSERT_FALSE(baseline_coords.empty());
  ASSERT_EQ(baseline.init_OBBTree(), MB_SUCCESS);

  const double scales[] = {1.0, 0.01, 2.5, 100.0};
  for (double scale : scales) {
    DagMC scaled;
    scaled.set_length_multiplier(scale);
    ASSERT_EQ(scaled.load_file(input_file), MB_SUCCESS);

    std::vector<double> scaled_coords = geometry_coordinates(scaled);
    ASSERT_EQ(scaled_coords.size(), baseline_coords.size());
    for (size_t i = 0; i < baseline_coords.size(); ++i)
      EXPECT_NEAR(scaled_coords[i], baseline_coords[i] * scale, 1e-10);
    EXPECT_DOUBLE_EQ(scaled.faceting_tolerance(),
                     baseline.faceting_tolerance() * scale);

    ASSERT_EQ(scaled.init_OBBTree(), MB_SUCCESS);
    EXPECT_NEAR(ray_distance(scaled), 5.0 * scale, 1e-8);
  }
}

TEST_F(DagmcSimpleTest, dagmc_length_multiplier_existing_contents) {
  const char scaled_input[] = "scaling_test.h5m";

  DagMC baseline;
  ASSERT_EQ(baseline.load_file(input_file), MB_SUCCESS);
  std::vector<double> baseline_coords = geometry_coordinates(baseline);
  ASSERT_EQ(baseline.init_OBBTree(), MB_SUCCESS);
  ASSERT_EQ(baseline.write_mesh(scaled_input, sizeof(scaled_input) - 1),
            MB_SUCCESS);

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  ASSERT_EQ(mbi->load_file(scaled_input), MB_SUCCESS);
  const double unrelated_coords[3] = {17.0, -23.0, 41.0};
  EntityHandle unrelated_vertex;
  ASSERT_EQ(mbi->create_vertex(unrelated_coords, unrelated_vertex), MB_SUCCESS);

  Tag resabs_tag;
  ASSERT_EQ(mbi->tag_get_handle("GEOMETRY_RESABS", 1, MB_TYPE_DOUBLE,
                                resabs_tag, MB_TAG_SPARSE | MB_TAG_CREAT),
            MB_SUCCESS);
  const EntityHandle root = 0;
  double root_resabs = 0.125;
  ASSERT_EQ(mbi->tag_set_data(resabs_tag, &root, 1, &root_resabs), MB_SUCCESS);
  EntityHandle unrelated_set;
  ASSERT_EQ(mbi->create_meshset(MESHSET_SET, unrelated_set), MB_SUCCESS);
  double unrelated_resabs = 9.0;
  ASSERT_EQ(mbi->tag_set_data(resabs_tag, &unrelated_set, 1, &unrelated_resabs),
            MB_SUCCESS);

  Tag obb_root_tag;
  ASSERT_EQ(mbi->tag_get_handle("OBB_ROOT", 1, MB_TYPE_HANDLE, obb_root_tag,
                                MB_TAG_SPARSE | MB_TAG_CREAT),
            MB_SUCCESS);
  EntityHandle unrelated_obb_set, unrelated_obb_root;
  ASSERT_EQ(mbi->create_meshset(MESHSET_SET, unrelated_obb_set), MB_SUCCESS);
  ASSERT_EQ(mbi->create_meshset(MESHSET_SET, unrelated_obb_root), MB_SUCCESS);
  ASSERT_EQ(mbi->tag_set_data(obb_root_tag, &unrelated_obb_set, 1,
                              &unrelated_obb_root),
            MB_SUCCESS);

  DagMC scaled(mbi);
  scaled.set_length_multiplier(0.4);
  ASSERT_EQ(scaled.load_existing_contents(), MB_SUCCESS);
  EXPECT_EQ(count_obb_roots(mbi.get()), 1u);
  EntityHandle current_obb_root;
  EXPECT_EQ(
      mbi->tag_get_data(obb_root_tag, &unrelated_obb_set, 1, &current_obb_root),
      MB_SUCCESS);
  EXPECT_EQ(current_obb_root, unrelated_obb_root);

  double current_root_resabs, current_unrelated_resabs;
  EXPECT_EQ(mbi->tag_get_data(resabs_tag, &root, 1, &current_root_resabs),
            MB_SUCCESS);
  EXPECT_EQ(mbi->tag_get_data(resabs_tag, &unrelated_set, 1,
                              &current_unrelated_resabs),
            MB_SUCCESS);
  EXPECT_DOUBLE_EQ(current_root_resabs, root_resabs * 0.4);
  EXPECT_DOUBLE_EQ(current_unrelated_resabs, unrelated_resabs);
  EXPECT_DOUBLE_EQ(scaled.faceting_tolerance(),
                   baseline.faceting_tolerance() * 0.4);

  std::vector<double> scaled_coords = geometry_coordinates(scaled);
  ASSERT_EQ(scaled_coords.size(), baseline_coords.size());
  for (size_t i = 0; i < baseline_coords.size(); ++i)
    EXPECT_NEAR(scaled_coords[i], baseline_coords[i] * 0.4, 1e-10);

  double current_unrelated_coords[3];
  ASSERT_EQ(mbi->get_coords(&unrelated_vertex, 1, current_unrelated_coords),
            MB_SUCCESS);
  for (int i = 0; i < 3; ++i)
    EXPECT_DOUBLE_EQ(current_unrelated_coords[i], unrelated_coords[i]);

  EXPECT_EQ(scaled.load_existing_contents(), MB_FAILURE);
  EXPECT_EQ(geometry_coordinates(scaled), scaled_coords);
  EXPECT_THROW(scaled.set_length_multiplier(2.0), std::logic_error);

  ASSERT_EQ(scaled.init_OBBTree(), MB_SUCCESS);
  EXPECT_NEAR(ray_distance(scaled), 2.0, 1e-8);
  remove(scaled_input);
}

TEST_F(DagmcSimpleTest, dagmc_length_multiplier_failed_load) {
  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  ASSERT_EQ(mbi->load_file(input_file), MB_SUCCESS);

  DagMC baseline(mbi);
  ASSERT_EQ(baseline.load_existing_contents(), MB_SUCCESS);
  std::vector<double> baseline_coords = geometry_coordinates(baseline);

  Range surfaces;
  ASSERT_EQ(baseline.geom_tool()->get_gsets_by_dimension(2, surfaces),
            MB_SUCCESS);
  ASSERT_FALSE(surfaces.empty());
  Tag obb_root_tag;
  ASSERT_EQ(mbi->tag_get_handle("OBB_ROOT", 1, MB_TYPE_HANDLE, obb_root_tag,
                                MB_TAG_SPARSE | MB_TAG_CREAT),
            MB_SUCCESS);
  Range vertices = geometry_vertices(baseline);
  ASSERT_FALSE(vertices.empty());
  EntityHandle invalid_root = vertices.front();
  ASSERT_EQ(
      mbi->tag_set_data(obb_root_tag, &surfaces.front(), 1, &invalid_root),
      MB_SUCCESS);

  DagMC scaled(mbi);
  scaled.set_length_multiplier(2.0);
  EXPECT_NE(scaled.load_existing_contents(), MB_SUCCESS);
  EXPECT_EQ(geometry_coordinates(baseline), baseline_coords);
  EXPECT_THROW(scaled.set_length_multiplier(3.0), std::logic_error);
  EXPECT_EQ(scaled.load_existing_contents(), MB_FAILURE);
}

TEST_F(DagmcSimpleTest, dagmc_length_multiplier_validation) {
  DagMC dagmc;
  EXPECT_DOUBLE_EQ(dagmc.length_multiplier(), 1.0);
  EXPECT_THROW(dagmc.set_length_multiplier(0.0), std::invalid_argument);
  EXPECT_THROW(dagmc.set_length_multiplier(-1.0), std::invalid_argument);
  EXPECT_THROW(
      dagmc.set_length_multiplier(std::numeric_limits<double>::quiet_NaN()),
      std::invalid_argument);
  EXPECT_THROW(
      dagmc.set_length_multiplier(std::numeric_limits<double>::infinity()),
      std::invalid_argument);
}

TEST_F(DagmcSimpleTest, dagmc_test_obb_retreval) {
  // make new dagmc
  std::cout << "test_obb_retreval" << std::endl;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();

  ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // write the file
  rval = dagmc->write_mesh("fcad", 4);

  dagmc.reset(new DagMC());
  rval = dagmc->load_file("fcad");
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // delete the fcad file
  remove("fcad");
}

TEST_F(DagmcSimpleTest, dagmc_build_obb) {
  ErrorCode rval = DAG->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_num_vols) {
  int expect_num_vols = 2;
  int num_vols = DAG->num_entities(3);
  EXPECT_EQ(expect_num_vols, num_vols);
}

TEST_F(DagmcSimpleTest, dagmc_point_in) {
  int result = 0;
  int expect_result = 1;
  int vol_idx = 1;
  double xyz[3] = {0.0, 0.0, 0.0};
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  ErrorCode rval = DAG->point_in_volume(vol_h, xyz, result);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_EQ(expect_result, result);
}

#if MOAB_VERSION_MAJOR == 5 && MOAB_VERSION_MINOR > 2
TEST_F(DagmcSimpleTest, dagmc_find_volume) {
  int vol_idx = 1;
  int vol_dim = 3;
  double xyz[3] = {0.0, 0.0, 0.0};
  EntityHandle vol_h;
  EntityHandle expected_vol_h = DAG->entity_by_index(vol_dim, vol_idx);
  ErrorCode rval = DAG->find_volume(xyz, vol_h);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_EQ(expected_vol_h, vol_h);
}
#endif

TEST_F(DagmcSimpleTest, dagmc_test_obb_retreval_rayfire) {
  // make new dagmc
  std::cout << "test_obb_retreval and ray_fire" << std::endl;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();

  ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // write the file
  rval = dagmc->write_mesh("fcad", 4);

  // now create new DAGMC
  dagmc.reset(new DagMC());
  rval = dagmc->load_file("fcad");
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // delete the fcad file
  remove("fcad");

  // now perform full ray fire
  double eps = 1.e-6;
  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  rval = DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_NEAR(expect_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcSimpleTest, dagmc_rayfire) {
  const double eps = 1e-6;  // epsilon for test, faceting tol?

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  ErrorCode rval = DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_NEAR(expect_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcSimpleTest, dagmc_closest_to) {
  const double eps = 1e-6;  // epsilon for test, faceting tolerance

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {-6.0, 0.0, 0.0};
  double distance;  // distance from point to nearest surface
  double expect_distance = 1.0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  ErrorCode rval = DAG->closest_to_location(vol_h, xyz, distance);
  EXPECT_EQ(rval, MB_SUCCESS);
  // distance should be 1.0 cm
  EXPECT_NEAR(expect_distance, distance, eps);
}

TEST_F(DagmcSimpleTest, dagmc_test_boundary) {
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  int surf_idx = 1;
  EntityHandle surf_h = DAG->entity_by_index(2, surf_idx);

  double xyz[3] = {0.0, 0.0, 5.0};
  double dir[3] = {0.0, 0.0, 1.0};
  int result;
  int expect_result = 0;

  ErrorCode rval = DAG->test_volume_boundary(vol_h, surf_h, xyz, dir, result);
  EXPECT_EQ(rval, MB_SUCCESS);
  // check ray leaving volume
  EXPECT_EQ(expect_result, result);
}

TEST_F(DagmcSimpleTest, dagmc_test_get_obb) {
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  double llc[3], urc[3];
  ErrorCode rval = DAG->getobb(vol_h, llc, urc);
  EXPECT_EQ(rval, MB_SUCCESS);

  // hardcoded value for 'test_geom.h5m'
  double geom_extent = 5.0;
  for (int i = 0; i < 3; i++) {
    EXPECT_LE(llc[i], -geom_extent);
    EXPECT_GE(urc[i], geom_extent);
  }
}
