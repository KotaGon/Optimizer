#ifndef _MODEL_H
#define _MODEL_H

namespace myOptimizeName
{
  class DataModelClass;
  class OptModelClass;
  class modelClass;

  class DataModelClass
  {
    private:

    public: 
      DataModelClass() = default;
      virtual ~DataModelClass(){}
      virtual void free()
      { }
  };

  class OptModelClass
  {
    private:

    public:  
      OptModelClass() = default;
      virtual ~OptModelClass(){}
  
      virtual void free()
      { }
  };

  class modelClass 
  {
    protected:
      DataModelClass *DataModel;
      OptModelClass *OptModel;
    public:
      modelClass() = default;
      virtual ~modelClass(){}
      
      virtual void importData(){};
      virtual void optimize(){};
      virtual void exportData(){};
      virtual void free()
      { 
	DataModel->free();
       	OptModel->free(); 
      
	delete DataModel;
	delete OptModel;
      }

      virtual DataModelClass *getDataModel(){ return DataModel; }
      virtual OptModelClass *getOptModel(){ return OptModel; }
  };
 

};

#endif
