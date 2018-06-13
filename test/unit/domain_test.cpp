/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include "fixtures/bestest-fixture.hpp"
#include "fixtures/typical-fixture.hpp"

#include "Domain.hpp"
#include "Errors.hpp"

using namespace Kiva;

static const double PI = 4.0*atan(1.0);


class DomainFixture : public BESTESTFixture {
protected:
  void SetUp() {
    specifySystem();
    ground = std::make_shared<Ground>(fnd,outputMap);
    fnd.createMeshData();
    domain = ground->domain;
    domain.setDomain(fnd);
  };

  Domain domain;
};

TEST_F( DomainFixture, domain_basics)
{
  EXPECT_EQ(domain.nX, 41);
  EXPECT_EQ(domain.nY, 1);
  EXPECT_EQ(domain.nZ, 19);
  EXPECT_EQ(domain.stepsize_i, 1);
  EXPECT_EQ(domain.stepsize_j, 41);
  EXPECT_EQ(domain.stepsize_k, 41);

  EXPECT_EQ(domain.dest_index_vector.size(), 3);
  EXPECT_EQ(domain.dest_index_vector[2].size(), domain.nX*domain.nY*domain.nZ);
}

TEST_F( DomainFixture, surface_indices)
{
  EXPECT_EQ(ground->foundation.surfaces[0].indices.size(), domain.nZ);
  EXPECT_EQ(ground->foundation.surfaces[4].indices.size(), domain.nX);
  EXPECT_EQ(ground->foundation.surfaces[5].indices.size(), 11);
}

TEST_F( DomainFixture, surface_tilt)
{
  EXPECT_DOUBLE_EQ(ground->foundation.surfaces[0].tilt, PI/2);
  EXPECT_DOUBLE_EQ(ground->foundation.surfaces[4].tilt, PI);
  EXPECT_DOUBLE_EQ(ground->foundation.surfaces[5].tilt, 0.0);
}

TEST_F( DomainFixture, cell_vector)
{
  EXPECT_EQ(domain.cell.size(), domain.nX*domain.nY*domain.nZ);

  EXPECT_EQ(domain.cell[0]->cellType, CellType::BOUNDARY);
  EXPECT_EQ(domain.cell[49]->cellType, CellType::NORMAL);
}

