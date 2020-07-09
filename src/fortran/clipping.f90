module sutherland_hodgman
    use iso_c_binding
    implicit none
    
    type point
        real(kind=8) x
        real(kind=8) y
    end type point
    
    type vector
        real(kind=8) x
        real(kind=8) y
    end type vector
    
    type intersection_result
        logical success
        type(point) p
    end type intersection_result

contains

    function cross(U, V) result(cp)
        type(vector), intent(in) :: U, V
        real(kind=8) :: cp
        cp = U%x * V%y - U%y * V%x
    end function cross
    

    function dot(U, V) result(dp)
        type(vector), intent(in) :: U, V
        real(kind=8) :: dp
        dp = U%x * V%x + U%y * V%y
    end function dot
   

    subroutine push(array, n, value)
        integer, intent(inout) :: n
        real(kind=8), dimension(:, :), intent(inout) :: array
        type(point), intent(in) :: value
        array(1, n + 1) = value%x
        array(2, n + 1) = value%y
        n = n + 1
    end subroutine push
    

    subroutine copyto(src, dst, n)
        real(kind=8), dimension(:, :), intent(in) :: src
        real(kind=8), dimension(:, :), intent(inout) :: dst
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            dst(1, i) = src(1, i)
            dst(2, i) = src(2, i)
        end do
    end subroutine copyto
    

    function inside(p, r, U) result(is_inside)
        type(point) :: p, r
        type(vector) :: U
        logical :: is_inside
        is_inside = U%x * (p%y - r%y) > U%y * (p%x - r%x)
    end function inside
    

    function intersection(a, V, r, N) result(int_res)
        type(intersection_result) :: int_res
        type(point) :: a, r, p
        type(vector) :: V, N, W
        real(kind=8) :: nw, nv, t
        W = vector(r%x - a%x, r%y - a%y)
        nw = dot(N, W)
        nv = dot(N, V)
        if (nv /= 0.0d0) then
            t = nw / nv
            p = point(a%x + t * V%x, a%y + t * V%y)
            int_res = intersection_result(.true., p)
        else
            p = point(0.0d0, 0.0d0)
            int_res = intersection_result(.false., p)
        end if   
    end function intersection
    

    function polygon_area(polygon, length) result(area)
        real(kind=8) :: area
        type(point) :: a, b, c
        type(vector) :: U, V
        real(kind=8), dimension(:, :) :: polygon
        integer :: length, i
        
        area = 0.0d0
        a = point(polygon(1, 1), polygon(2, 1))
        b = point(polygon(1, 2), polygon(2, 2))
        U = vector(b%x - a%x, b%y - a%y)
        do i = 3, length
            c = point(polygon(1, i), polygon(2, i))
            V = vector(a%x - c%x, a%y - c%y)
            area = area + abs(cross(U, V))
            b = c
            U = V
        end do
        area = 0.5d0 * area
    end function polygon_area
    

    function clip_polygons(polygon, clipper, n_max) result(area)
        type(point) :: r, s, a, b, p
        type(vector) :: U, N, V
        type(intersection_result) :: int_res
        integer, dimension(2) :: polygon_shape, clipper_shape
        integer :: n_output, n_poly, n_clip, length, i, j
        integer(kind=8) :: n_max
        real(kind=8), dimension(:, :) :: polygon, clipper
        real(kind=8), dimension(2, n_max):: output, subject
        real(kind=8) :: area
        logical :: a_inside, b_inside
        
        polygon_shape = shape(clipper)
        clipper_shape = shape(polygon)
        n_poly = polygon_shape(2)
        n_clip = clipper_shape(2)

        n_output = n_poly
        call copyto(polygon, output, n_output)

        r = point(clipper(1, n_clip), clipper(2, n_clip))
        do i = 1, n_clip
            
            s = point(clipper(1, i), clipper(2, i))
            U = vector(s%x - r%x, s%y - r%y)
            N = vector(-U%y, U%x)
            
            if (U%x == 0.0d0 .and. U%y == 0.0d0) then
                cycle
            end if
            
            length = n_output
            call copyto(output, subject, length)
            n_output = 0
            
            a = point(subject(1, length), subject(2, length))
            a_inside = inside(a, r, U)
            do j = 1, length
                b = point(subject(1, j), subject(2, j))
                V = vector(b%x - a%x, b%y - a%y)
                
                if (V%x == 0.0d0 .and. V%y == 0.0d0) then
                    cycle
                end if
                
                b_inside = inside(b, r, U)
                
                if (b_inside) then
                    if (.not. a_inside) then
                        int_res = intersection(a, V, r, N)
                        if (int_res%success) then
                            call push(output, n_output, int_res%p)
                        end if
                    end if
                    call push(output, n_output, b)
                else if (a_inside) then
                    int_res = intersection(a, V, r, N)
                    if (int_res%success) then
                        call push(output, n_output, int_res%p)
                    else
                        b_inside = .true.
                        call push(output, n_output, b)
                    end if
                end if
                a = b
                a_inside = b_inside
            end do
            
            if (n_output < 3) then
                area = 0.0d0
                return
            end if
            
            r = s
        end do
        area = polygon_area(output, n_output)
    end function clip_polygons
    

    subroutine area_of_intersection(ndim, nvertex, ntriangles, polygons, clippers, areas) bind(c, name="area_of_intersection")
        !DEC$ ATTRIBUTES DLLEXPORT :: area_of_intersection
        integer(c_int64_t) :: ntriangles, nvertex, ndim, i
        real(c_double), dimension(ndim, nvertex, ntriangles), intent(in) :: polygons, clippers
        real(c_double), dimension(ntriangles), intent(out) :: areas
        do i = 1, ntriangles
            areas(i) = clip_polygons(polygons(:, :, i), clippers(:, :, i), nvertex * 2)
        end do
    end subroutine area_of_intersection

end