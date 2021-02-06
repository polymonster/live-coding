void basis_from_axis(const vec3d axis, vec3d& right, vec3d& up, vec3d& at)
{
    right = cross(axis, vec3d::unit_y());
    
    if (mag(right) < 0.1)
        right = cross(axis, vec3d::unit_z());
    
    if (mag(right) < 0.1)
        right = cross(axis, vec3d::unit_x());
        
    normalise(right);
    up = normalised(cross(axis, right));
    right = normalised(cross(axis, up));
    at = cross(right, up);
}

void pentagon_in_axis(const vec3d axis, const vec3d pos, f64 start_angle, bool recurse)
{
    vec3d right, up, at;
    basis_from_axis(axis, right, up, at);
    
    f64 half_gr = 1.61803398875l/2.0;
        
    f64 internal_angle = 0.309017 * 1.5;
    f64 angle_step = M_PI / 2.5;
    f64 a = start_angle;
    for(u32 i = 0; i < 5; ++i)
    {
        f64 x = sin(a) * M_INV_PHI;
        f64 y = cos(a) * M_INV_PHI;
        
        vec3d p = pos + right * x + up * y;
        
        a += angle_step;
        f64 x2 = sin(a) * M_INV_PHI;
        f64 y2 = cos(a) * M_INV_PHI;
        
        vec3d np = pos + right * x2 + up * y2;
        add_line((vec3f)p, (vec3f)np, vec4f::green());
                    
        vec3d ev = normalised(np - p);
        vec3d cp = normalised(cross(ev, axis));

        vec3d mid = p + (np - p) * 0.5;
        
        f64 rx = sin((M_PI*2.0)+internal_angle) * M_INV_PHI;
        f64 ry = cos((M_PI*2.0)+internal_angle) * M_INV_PHI;
        vec3d xp = mid + cp * rx + axis * ry;
        
        vec3d xv = normalised(xp - mid);

        if(recurse)
        {
            vec3d next_axis = normalised(cross(xv, ev));
            pentagon_in_axis(next_axis, mid + xv * half_gr * M_INV_PHI, M_PI + start_angle, false);
        }
    }
}

void penatgon_icosa(const vec3d axis, const vec3d pos, f64 start_angle)
{
    vec3d right, up, at;
    basis_from_axis(axis, right, up, at);
    
    vec3d tip = pos - at * M_INV_PHI;
    vec3d dip = pos + at * 0.5 * 2.0;
    
    f64 angle_step = M_PI / 2.5;
    f64 a = start_angle;
    for(u32 i = 0; i < 5; ++i)
    {
        f64 x = sin(a);
        f64 y = cos(a);
        
        vec3d p = pos + right * x + up * y;
        
        a += angle_step;
        f64 x2 = sin(a);
        f64 y2 = cos(a);
        
        vec3d np = pos + right * x2 + up * y2;
        add_line((vec3f)p, (vec3f)np, vec4f::magenta());
        add_line((vec3f)p, (vec3f)tip, vec4f::magenta());
        add_line((vec3f)np, (vec3f)tip, vec4f::magenta());
        
        vec3d side_dip = dip + cross(normalized(p-np), at);
        add_line((vec3f)np, (vec3f)side_dip, vec4f::magenta());
        add_line((vec3f)p, (vec3f)side_dip, vec4f::magenta());
    }
}

void icosahedron(vec3f axis, vec3f pos)
{
    penatgon_icosa((vec3d)axis, (vec3d)(pos + axis * 0.5f), 0.0);
    penatgon_icosa((vec3d)-axis, (vec3d)(pos - axis * 0.5f), M_PI);
}

void dodecahedron(vec3f axis, vec3f pos)
{
    f32 h = M_PI*0.83333333333f * 0.5f * M_INV_PHI;
    pentagon_in_axis((vec3d)axis, (vec3d)pos + vec3d(0.0, -h, 0.0), 0.0f, true);
    pentagon_in_axis((vec3d)-axis, (vec3d)pos + vec3d(0.0, h, 0.0), M_PI, true);
}

void tetrahedron(vec3d axis, vec3d pos)
{
    vertex_model* vertices = nullptr;
    
    vec3d right, up, at;
    basis_from_axis(axis, right, up, at);
        
    vec3d tip = pos - at * sqrt(2.0); // sqrt 2 is pythagoras constant
    
    f64 angle_step = (M_PI*2.0) / 3.0;
    f64 a = 0.0f;
    for(u32 i = 0; i < 3; ++i)
    {
        f64 x = sin(a);
        f64 y = cos(a);
                    
        vec3d p = pos + right * x + up * y;
        
        a += angle_step;
        f64 x2 = sin(a);
        f64 y2 = cos(a);
        
        vec3d np = pos + right * x2 + up * y2;
        
        vec3f n = maths::get_normal((vec3f)p, (vec3f)np, (vec3f)tip);
        vec3f b = (vec3f)normalised(p - np);
        vec3f t = cross(n, b);
        
        vertex_model v[3];
        
        v[0].pos.xyz = (vec3f)p;
        v[1].pos.xyz = (vec3f)np;
        v[2].pos.xyz = (vec3f)tip;
        
        for(u32 j = 0; j < 3; ++j)
        {
            v[j].pos.w = 1.0;
            v[j].normal = vec4f(n, 1.0f);
            v[j].bitangent = vec4f(b, 1.0f);
            v[j].tangent = vec4f(t, 1.0f);
        }

        add_triangle(v[0].pos.xyz, v[1].pos.xyz, v[2].pos.xyz, vec4f::cyan());
    }
    
    u16* indices = nullptr;
    
    vec3f min_extents = vec3f::flt_max();
    vec3f max_extents = -vec3f::flt_max();
    
    u32 nv = sb_count(vertices);
    for(u32 i = 0; i < nv; i++)
    {
        sb_push(indices, i);
        
        min_extents = min_union(vertices[i].pos.xyz, min_extents);
        max_extents = max_union(vertices[i].pos.xyz, max_extents);
    }

    sb_free(vertices);
}

void octahedron()
{
    vertex_model* vertices = nullptr;
    
    vec3f corner[] = {
        vec3f(-1.0, 0.0, -1.0),
        vec3f(-1.0, 0.0, 1.0),
        vec3f(1.0, 0.0, 1.0),
        vec3f(1.0, 0.0, -1.0)
    };
    
    f32 pc = sqrt(2.0);
    vec3f tip = vec3f(0.0f, pc, 0.0f);
    vec3f dip = vec3f(0.0f, -pc, 0.0f);
    
    for(u32 i = 0; i < 4; ++i)
    {
        u32 n = (i + 1) % 4;
        
        vec3f y[] = {
            tip,
            dip
        };
        
        // 2 tris per edg
        for(u32 j = 0; j < 2; ++j)
        {
            vertex_model v[3];
            v[0].pos.xyz = corner[i];

            v[1].pos.xyz = corner[n];
            v[2].pos.xyz = y[j];
            
            vec3f n = maths::get_normal(v[0].pos.xyz, v[1].pos.xyz, v[2].pos.xyz);
            vec3f b = normalised(v[0].pos.xyz - v[1].pos.xyz);
            vec3f t = cross(n, b);
            
            for(u32 k = 0; k < 3; ++k)
            {
                v[k].pos.w = 1.0f;
                v[k].normal = vec4f(n, 1.0f);
                v[k].tangent = vec4f(t, 1.0f);
                v[k].bitangent = vec4f(b, 1.0f);
                
                sb_push(vertices, v[k]);
            }
        }
    }
    
    u16* indices = nullptr;
    
    vec3f min_extents = vec3f::flt_max();
    vec3f max_extents = -vec3f::flt_max();
    
    u32 nv = sb_count(vertices);
    for(u32 i = 0; i < nv; i++)
    {
        sb_push(indices, i);
        
        min_extents = min_union(vertices[i].pos.xyz, min_extents);
        max_extents = max_union(vertices[i].pos.xyz, max_extents);
    }
    
    for(u32 i = 0; i < nv; i+=3)
    {
        dbg::add_triangle(vertices[i].pos.xyz, vertices[i+1].pos.xyz, vertices[i+2].pos.xyz);
    }
}

void create_torus_primitive(f32 radius)
{
    f64 angle_step = (M_PI*2.0)/64.0;
    f64 aa = 0.0f;
    for(u32 i = 0; i < 64; ++i)
    {
        f64 x = sin(aa);
        f64 y = cos(aa);
        
        aa += angle_step;
        f64 x2 = sin(aa);
        f64 y2 = cos(aa);
        
        f64 x3 = sin(aa + angle_step);
        f64 y3 = cos(aa + angle_step);
        
        vec3f p = vec3f(x, 0.0, y);
        vec3f np = vec3f(x2, 0.0, y2);
        vec3f nnp = vec3f(x3, 0.0, y3);
        
        vec3f at = normalized(np - p);
        vec3f up = vec3f::unit_y();
        vec3f right = cross(up, at);
        
        vec3f nat = normalized(nnp - np);
        vec3f nright = cross(up, nat);
        
        f64 ab = 0.0f;
        for(u32 j = 0; j < 64; ++j)
        {
            f32 vx = sin(ab) * radius;
            f32 vy = cos(ab) * radius;

            vec3f vv = p + vx * up + vy * right;          
            ab += angle_step;         

            f32 vx2 = sin(ab) * radius;
            f32 vy2 = cos(ab) * radius;
            
            vec3f vv2 = p + vx2 * up + vy2 * right;
            
            add_line(vv, vv2);
            
            vec3f vv3 = np + vx * up + vy * nright;
            vec3f vv4 = np + vx2 * up + vy2 * nright;
            
            //add_line(vv, vv2);
            
            //add_line(vv, vv3, vec4f::yellow());
            
            add_line(vv, vv3, vec4f::yellow());
        }
    }
}

void prims()
{
    tetrahedron(vec3d::unit_y(), vec3d(-2.5f, 0.0f, 0.0f));
    octahedron();
    dodecahedron(vec3f::unit_y(), vec3f(2.5f, 0.0f, 0.0f));
    icosahedron(vec3f::unit_y(), vec3f(5.0f, 0.0f, 0.0f));
}