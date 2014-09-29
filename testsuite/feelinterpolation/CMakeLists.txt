
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
set_directory_properties(PROPERTIES LABEL testinterpolation )

feelpp_add_test( form_interpolation NO_MPI_TEST TIMEOUT 40)
feelpp_add_test( interpolation NO_MPI_TEST   )
feelpp_add_test( interpolation_op NO_MPI_TEST TIMEOUT 60 )
feelpp_add_test( interp_twomesh NO_MPI_TEST )
feelpp_add_test( eval_at_point TIMEOUT 60)
feelpp_add_test( operatorinterpolation )
feelpp_add_test( hypercubeinterpolation )

feelpp_add_test( interpolation_RT CFG test_interpolation_RT.cfg NO_MPI_TEST)
feelpp_add_test( interpolation_Nedelec CFG test_interpolation_Nedelec.cfg NO_MPI_TEST)
feelpp_add_test( interpolation_Nedelec3d CFG test_interpolation_Nedelec3d.cfg NO_MPI_TEST TIMEOUT 40)
