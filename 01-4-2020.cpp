#include "cr/cr.h"
#include "ecs/ecs_scene.h"
#include "ecs/ecs_utilities.h"
#include "ecs/ecs_resources.h"

#define DLL 1
#include "ecs/ecs_live.h"
#include "str/Str.cpp"

#include "renderer.h"
#include "maths/maths.h"
#include "../../shader_structs/forward_render.h"

#include <stdio.h>

using namespace pen;
using namespace put;
using namespace ecs;
using namespace dbg;

s32 randri(s32 min, s32 max)
{
    s32 range = max - min;
    if(range == 0)
        return min;
    return min + rand()%range;
}

struct live_lib : public __ecs, public __dbg
{
    const f32           cube_pad = 0.5f;
    const f32           cube_size = 10.0f;
    const f32           grid_cube_size = cube_size*2;
    const vec3f         base_offset = vec3f(-(f32)grid_size*cube_size, 0.0f, -(f32)grid_size*cube_size);
    static const u32    grid_size = 50;
    static const u32    num_boxes = grid_size*grid_size;
    
    ecs_scene*          scene;
    u32                 occ[grid_size][grid_size];
    u32                 box_start;
    
    u32 box_count = 0;
    
    void init(live_context* ctx)
    {
        scene = ctx->scene;
        memcpy(&__ecs_start, &ctx->ecs_funcs->__ecs_start, (intptr_t)&ctx->ecs_funcs->__ecs_end - (intptr_t)&ctx->ecs_funcs->__ecs_start);
        memcpy(&__dbg_start, &ctx->dbg_funcs->__dbg_start, (intptr_t)&ctx->dbg_funcs->__dbg_end - (intptr_t)&ctx->dbg_funcs->__dbg_start);
    }
    
    int on_load(live_context* ctx)
    {
        init(ctx);
        
        material_resource* default_material = get_material_resource(PEN_HASH("default_material"));
        geometry_resource* box_resource = get_geometry_resource(PEN_HASH("cube"));
            
        clear_scene(scene);

        // directional light
        u32 light = get_new_entity(scene);
        instantiate_light(scene, light);
        scene->names[light] = "front_light";
        scene->id_name[light] = PEN_HASH("front_light");
        scene->lights[light].colour = vec3f(0.3f, 0.5f, 0.8f);
        scene->lights[light].direction = normalised(vec3f(-1.0f, 0.5f, -0.5f));
        scene->lights[light].type = e_light_type::dir;
        scene->lights[light].flags |= e_light_flags::shadow_map;
        scene->transforms[light].translation = vec3f::zero();
        scene->transforms[light].rotation = quat();
        scene->transforms[light].scale = vec3f::one();
        scene->entities[light] |= e_cmp::light;
        scene->entities[light] |= e_cmp::transform;
        
        light = get_new_entity(scene);
        instantiate_light(scene, light);
        scene->names[light] = "front_light";
        scene->id_name[light] = PEN_HASH("front_light");
        scene->lights[light].colour = vec3f(1.0f, 0.5f, 0.3f);
        scene->lights[light].direction = normalised(vec3f(-1.0f, 0.3f, -0.8f));
        scene->lights[light].type = e_light_type::dir;
        scene->lights[light].flags |= e_light_flags::shadow_map;
        scene->transforms[light].translation = vec3f::zero();
        scene->transforms[light].rotation = quat();
        scene->transforms[light].scale = vec3f::one();
        scene->entities[light] |= e_cmp::light;
        scene->entities[light] |= e_cmp::transform;
        
        light = get_new_entity(scene);
        instantiate_light(scene, light);
        scene->names[light] = "front_light";
        scene->id_name[light] = PEN_HASH("front_light");
        scene->lights[light].colour = vec3f(0.3f, 0.4f, 0.4f);
        scene->lights[light].direction = normalised(vec3f(-1.0f, 0.8f, 0.8f));
        scene->lights[light].type = e_light_type::dir;
        scene->lights[light].flags |= e_light_flags::shadow_map;
        scene->transforms[light].translation = vec3f::zero();
        scene->transforms[light].rotation = quat();
        scene->transforms[light].scale = vec3f::one();
        scene->entities[light] |= e_cmp::light;
        scene->entities[light] |= e_cmp::transform;
        
        // master instance
        u32 master = get_new_entity(scene);
        scene->names[master] = "master";
        scene->id_name[master] = PEN_HASH("master");
        scene->transforms[master].rotation = quat();
        scene->transforms[master].scale = vec3f::one();
        scene->transforms[master].translation = vec3f::zero();
               
        memset(&occ[0][0], 0x00, sizeof(u32)*grid_size*grid_size);
        srand(2000);
        u32 cc = 1;
        for(;;)
        {
            // find first free block
            vec2i start;
            for(u32 i = 0; i < grid_size; ++i)
                for(u32 j = 0; j < grid_size; ++j)
                    if(occ[i][j] == 0)
                    {
                        start = vec2i(j, i);
                        goto first_block;
                    }
            break;
            first_block:
            
            // find max range
            vec2i range_min = vec2i(2, 2);
            vec2i range = vec2i(5, 5);
            for(u32 i = start.x; i < start.x+range.x; ++i)
                if(i == grid_size || occ[start.y][i] != 0)
                {
                    range.x = i-start.x;
                    break;
                }

            for(u32 i = start.y; i < start.y+range.y; ++i)
                if(i == grid_size || occ[i][start.x] != 0)
                {
                    range.y = i-start.y-1;
                    break;
                }
                                      
            range_min = max_union(min_union(range_min, range), vec2i(1));
            
            
            // randomise size x and y, get pos and scale of cubes
            vec2i size = vec2i(randri(range_min.x, range.x), randri(range_min.y, range.y));
            vec3f start_pos = vec3f(start.x, 0.0f, start.y) * cube_size*2.0f;
            vec3f scale = vec3f(size.x, 1.0f, size.y) * cube_size;
            vec3f pos = start_pos + scale;
            
            if(range_min.x > range.x || range_min.x > range.x)
                pos.y += 1.0f + cube_size * 0.75f * ((f32)box_count/num_boxes);
            
            vec3f hsv = vec3f(lerp(0.0f, 0.2f, (f32)start.x/grid_size) + 0.2f, 1.0f, 0.5f);
                        
            u32 box = get_new_entity(scene);
            scene->names[box] = "box";
            scene->id_name[box] = PEN_HASH("box");
            scene->transforms[box].rotation = quat();
            scene->transforms[box].scale = scale - cube_pad;
            scene->transforms[box].translation = base_offset + pos;
            scene->parents[box] = box;
            scene->entities[box] |= e_cmp::transform;
            instantiate_geometry(box_resource, scene, box);
            instantiate_material(default_material, scene, box);
            instantiate_model_cbuffer(scene, box);
            
            forward_render::forward_lit* mat = (forward_render::forward_lit*)&scene->material_data[box].data[0];
            mat->m_albedo = vec4f(hsv, 1.0f);
                
            if(box_count == 0)
                box_start = box;
            box_count++;
            
            // mark occupied
            for(u32 i = start.y; i < start.y + size.y; ++i)
                for(u32 j = start.x; j < start.x + size.x; ++j)
                    occ[i][j] = cc;
        }
        
        instance_entity_range(scene, master, box_count);
   
        return 0;
    }
    
    void debug_grid()
    {
        for(u32 i = 0; i < grid_size; ++i)
            for(u32 j = 0; j < grid_size; ++j)
            {
                vec4f col = vec4f::green();
                if(occ[i][j] == 0)
                    col = vec4f::red();
                
                vec3f bbmin = base_offset + vec3f(j*grid_cube_size, 0.0f, i*grid_cube_size);
                add_aabb(bbmin, bbmin + vec3f(grid_cube_size), col);
            }
    }
    
    int on_update(f32 dt)
    {
        static f32 p = 0.0f;
        static vec2f wp = base_offset.xz;
        
        for(u32 i = 0; i < box_count; ++i)
        {
            f32 dz = floor(min(scene->transforms[i].translation.x - wp.x, 0.0f)) / 300.0f;
            f32 dx = floor(min(scene->transforms[i].translation.z - wp.y, 0.0f)) / 300.0f;
            
            f32 s = scene->transforms[i].scale.x * 0.5f;
            f32 w = sin(-dz-dx) * s;
            
            scene->transforms[i].scale.y = max(w, 10.0f);
            scene->entities[i] |= e_cmp::transform;
        }
        
        wp += dt * 100.0f;

        p += dt;
        return 0;
    }
    
    int on_unload()
    {
        return 0;
    }
};

CR_EXPORT int cr_main(struct cr_plugin *ctx, enum cr_op operation)
{
    live_context* live_ctx = (live_context*)ctx->userdata;
    static live_lib ll;
    
    switch (operation)
    {
        case CR_LOAD:
            return ll.on_load(live_ctx);
        case CR_UNLOAD:
            return ll.on_unload();
        case CR_CLOSE:
            return 0;
        default:
            break;
    }
    
    return ll.on_update(live_ctx->dt);
}

namespace pen
{
    pen_creation_params pen_entry(int argc, char** argv)
    {
        return {};
    }
    
    void* user_entry(void* params)
    {
        return nullptr;
    }
}


