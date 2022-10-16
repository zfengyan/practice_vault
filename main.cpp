#include <iostream>
#include <functional> // for std::reference_wrapper<T>
#include <future> // for std::async
#include <mutex> // for std::mutex
#include <vector>

struct Mesh;
class LoadMesh;

template<typename T>
using Ref = std::reference_wrapper<T>;

struct Mesh {
    static Ref<Mesh> Load(const std::string& file) {
        // ...
        Mesh m;
        return std::ref(m); // can not use std::ref(Mesh())

        // generate an object of type std::reference_wrapper
    }
};


namespace GenerateMesh {

    std::vector<Ref<Mesh>> m_Meshes;
    std::vector<std::future<void>> m_Futures; // store the return value of std::async, necessary step to make async work
    std::mutex s_MeshesMutex;


    void LoadMesh(std::vector<Ref<Mesh>>* meshes, std::string filepath /* make copy */) {
        auto mesh = Mesh::Load(filepath);

        // using a local lock_guard to lock mtx guarantees unlocking on destruction / exception:
        std::lock_guard<std::mutex> lock(s_MeshesMutex); // lock the meshes to avoid conflict
        meshes->emplace_back(mesh);

        //std::cout<<mesh;
    }


    void LoadMeshes() {
        std::vector<std::string> meshFilepaths; // read from a .txt file ...

        // ...
        // meshFilepaths.emplace_back(...)

        /*for(const auto& file : meshFilepaths){
            m_Meshes.emplace_back(Mesh::Load(file));
        }*/

        for (const auto& file : meshFilepaths) {
            m_Futures.emplace_back(std::async(std::launch::async, LoadMesh, &m_Meshes, file)); /* using auto result = std::async(...) ? */
            // question: do we have to use get()?
            // i.e. 
            // auto result = std::async(std::launch::async, LoadMesh, &m_Meshes, file);
            // result.get();
            // ?
        }
    }

};


int main()
{
    GenerateMesh::LoadMeshes();
    std::cout << "done" << '\n';
    return 0;
}