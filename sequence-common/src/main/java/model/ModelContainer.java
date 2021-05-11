package model;

import java.util.*;

public class ModelContainer implements Iterable<ModelInterface> {
    private final Map<Class, ModelInterface> map = new HashMap<>();

    public void addAll(List<ModelInterface> models) {
        models.forEach(this::add);
    }

    @SuppressWarnings("unchecked")
    public <T> T get(Class<T> klass) {
        return (T) map.get(klass);
    }

    public boolean contains(Class klass) {
        return map.containsKey(klass);

    }

    public void add(ModelInterface model) {
        this.map.put(model.getClass(), model);

    }

    @Override
    public Iterator<ModelInterface> iterator() {
        return map.values().iterator();
    }
}
